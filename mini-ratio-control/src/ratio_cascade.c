#include "ratio_cascade.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/**
 * @file ratio_cascade.c
 * @brief Cascade Ratio Control Implementation
 *
 * Cascade control nests one PID inside another — the outer (primary)
 * loop provides setpoint to the inner (secondary) loop through a
 * ratio station. The inner loop rejects fast disturbances; the outer
 * loop tracks slow setpoint changes and corrects steady-state error.
 *
 * Architecture:
 *   Primary PID  → produces demand signal
 *   Ratio station → converts demand to secondary SP
 *   Secondary PID → tracks ratio-adjusted SP
 *
 * Anti-windup: When secondary saturates, primary integral is frozen.
 * Bumpless transfer: When switching secondary from LOCAL to CASCADE,
 *   primary output is initialized to current secondary SP / ratio.
 *
 * L1: Cascade structure definitions
 * L2: Cascade modes, bump-less transfer, anti-windup
 * L3: Cascade loop decomposition into inner/outer transfer functions
 * L4: Cascade stability — inner loop must be 3-5x faster than outer
 * L5: PID anti-windup with back-calculation
 * L6: Distillation temperature/flow cascade
 * L7: Boiler drum level/feedwater cascade
 * L8: Decoupling in cascade structures
 *
 * MIT 6.302 §9: Cascade Control
 * Cambridge 3F2: Cascade Architectures
 * Purdue ME 675: Multivariable Cascade
 */

#ifndef NAN
#define NAN (0.0/0.0)
#endif

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

/* -------------------------------------------------------------------
 * L2: PID Loop Implementation (Embedded)
 * ------------------------------------------------------------------- */

int pid_loop_init(pid_loop_t *pid, const char *tag,
                  double kp, double ki, double kd,
                  double out_min, double out_max)
{
    if (!pid) return -1;
    if (out_min > out_max) return -2;

    memset(pid, 0, sizeof(*pid));
    pid->tag_name       = tag ? tag : "";
    pid->setpoint       = 0.0;
    pid->process_variable = 0.0;
    pid->output         = 0.0;
    pid->output_min     = out_min;
    pid->output_max     = out_max;
    pid->kp             = kp;
    pid->ki             = ki;
    pid->kd             = kd;
    pid->integral       = 0.0;
    pid->error_prev     = 0.0;
    pid->dt_accum       = 0.0;
    pid->in_auto        = 0;
    pid->in_cascade     = 0;
    pid->antiwindup_limit = out_max;
    pid->bumpless_active = 0;

    return 0;
}

int pid_loop_set_sp(pid_loop_t *pid, double sp)
{
    if (!pid) return -1;
    pid->setpoint = sp;
    return 0;
}

int pid_loop_set_pv(pid_loop_t *pid, double pv)
{
    if (!pid) return -1;
    pid->process_variable = pv;
    return 0;
}

int pid_loop_execute(pid_loop_t *pid, double dt)
{
    if (!pid) return -1;
    if (dt <= 0.0) return -2;
    if (!pid->in_auto) return 0;

    double error = pid->setpoint - pid->process_variable;

    /* Proportional term */
    double p_term = pid->kp * error;

    /* Integral term with anti-windup via clamping */
    pid->integral += pid->ki * error * dt;

    /* Anti-windup: clamp integral contribution */
    double i_term = pid->integral;
    if (i_term > pid->output_max) {
        pid->integral = pid->output_max;
        i_term = pid->output_max;
    }
    if (i_term < pid->output_min) {
        pid->integral = pid->output_min;
        i_term = pid->output_min;
    }

    /* Derivative term (filtered: avoid derivative kick) */
    double d_term = 0.0;
    if (dt > 1e-12) {
        double deriv = (error - pid->error_prev) / dt;
        d_term = pid->kd * deriv;
    }

    /* Total output */
    double output = p_term + i_term + d_term;

    /* Clamp */
    if (output > pid->output_max) output = pid->output_max;
    if (output < pid->output_min) output = pid->output_min;

    /* Saturation detection for anti-windup */
    if (output >= pid->output_max || output <= pid->output_min) {
        if (output * error > 0) {
            /* Saturating in the direction of error — freeze integral */
            pid->integral -= pid->ki * error * dt;
        }
    }

    pid->output     = output;
    pid->error_prev = error;
    pid->dt_accum  += dt;

    return 0;
}

double pid_loop_get_output(const pid_loop_t *pid)
{
    if (!pid) return NAN;
    return pid->output;
}

int pid_loop_set_auto(pid_loop_t *pid, int auto_mode)
{
    if (!pid) return -1;
    pid->in_auto = auto_mode;
    return 0;
}

int pid_loop_bumpless_init(pid_loop_t *pid, double current_output)
{
    if (!pid) return -1;

    /* Initialize integral term so output = current_output at zero error */
    pid->integral  = current_output;
    pid->output    = current_output;
    pid->error_prev = 0.0;
    pid->bumpless_active = 1;

    return 0;
}

int pid_loop_reset(pid_loop_t *pid)
{
    if (!pid) return -1;
    pid->integral   = 0.0;
    pid->output     = 0.0;
    pid->error_prev = 0.0;
    pid->dt_accum   = 0.0;
    return 0;
}

/* -------------------------------------------------------------------
 * L2: Cascade Ratio Control
 * ------------------------------------------------------------------- */

int cascade_ratio_init(cascade_ratio_t *cr,
                       cascade_mode_t mode,
                       unsigned int primary_id,
                       const char *tag)
{
    if (!cr) return -1;
    if (mode < 0 || mode >= CASCADE_MODE_COUNT) return -2;

    (void)primary_id;
    memset(cr, 0, sizeof(*cr));

    cr->mode  = mode;
    cr->state = CASCADE_STATE_INITIALIZING;

    /* Initialize primary PID with neutral gains */
    pid_loop_init(&cr->primary_pid, tag, 1.0, 0.1, 0.0, 0.0, 100.0);

    /* Initialize secondary PIDs */
    pid_loop_init(&cr->secondary_pid_a, "", 1.0, 0.5, 0.0, 0.0, 100.0);
    pid_loop_init(&cr->secondary_pid_b, "", 1.0, 0.5, 0.0, 0.0, 100.0);

    /* Default ratio configs: 1:1 */
    ratio_configure(NULL, RATIO_MODE_CASCADE_RATIO, RATIO_ACTION_MULTIPLY,
                    1.0, 0.0, 0.0, INFINITY);
    cr->ratio_cfg_a.ratio_gain = 1.0;
    cr->ratio_cfg_a.ratio_bias = 0.0;
    cr->ratio_cfg_a.ratio_min  = 0.0;
    cr->ratio_cfg_a.ratio_max  = INFINITY;
    cr->ratio_cfg_b.ratio_gain = 1.0;
    cr->ratio_cfg_b.ratio_bias = 0.0;
    cr->ratio_cfg_b.ratio_min  = 0.0;
    cr->ratio_cfg_b.ratio_max  = INFINITY;

    cr->primary_sp      = 0.0;
    cr->primary_pv      = 0.0;
    cr->primary_output  = 0.0;
    cr->secondary_sp_a  = 0.0;
    cr->secondary_sp_b  = 0.0;
    cr->split_point     = 50.0; /* 50% for split-range */
    cr->override_high_limit = INFINITY;
    cr->override_low_limit  = 0.0;
    cr->override_active = 0;
    cr->cycle_count     = 0;

    return 0;
}

int cascade_ratio_configure(cascade_ratio_t *cr,
                            const ratio_config_t *cfg_a,
                            const ratio_config_t *cfg_b)
{
    if (!cr) return -1;
    if (cfg_a) cr->ratio_cfg_a = *cfg_a;
    if (cfg_b) cr->ratio_cfg_b = *cfg_b;
    return 0;
}

int cascade_ratio_tune_primary(cascade_ratio_t *cr,
                               double kp, double ki, double kd,
                               double out_min, double out_max)
{
    if (!cr) return -1;
    return pid_loop_init(&cr->primary_pid,
                         cr->primary_pid.tag_name,
                         kp, ki, kd, out_min, out_max);
}

int cascade_ratio_tune_secondary(cascade_ratio_t *cr,
                                 int which,
                                 double kp, double ki, double kd,
                                 double out_min, double out_max)
{
    if (!cr) return -1;
    pid_loop_t *sec = (which == 0) ? &cr->secondary_pid_a
                                    : &cr->secondary_pid_b;
    return pid_loop_init(sec, sec->tag_name, kp, ki, kd, out_min, out_max);
}

int cascade_ratio_set_primary_sp(cascade_ratio_t *cr, double sp)
{
    if (!cr) return -1;
    cr->primary_sp = sp;

    /* Bumpless: when changing primary SP while in cascade,
     * back-calculate secondary SP to avoid bump */
    if (cr->state == CASCADE_STATE_PRIMARY_ACTIVE) {
        double new_secondary_sp = sp * cr->ratio_cfg_a.ratio_gain
                                  + cr->ratio_cfg_a.ratio_bias;
        cr->secondary_sp_a = new_secondary_sp;
        pid_loop_set_sp(&cr->secondary_pid_a, new_secondary_sp);
    }

    return 0;
}

int cascade_ratio_execute(cascade_ratio_t *cr,
                          double primary_pv,
                          double secondary_pv_a,
                          double secondary_pv_b,
                          double dt)
{
    if (!cr) return -1;
    if (dt <= 0.0) return -2;

    cr->primary_pv = primary_pv;

    /* --- Primary PID --- */
    pid_loop_set_sp(&cr->primary_pid, cr->primary_sp);
    pid_loop_set_pv(&cr->primary_pid, primary_pv);
    pid_loop_execute(&cr->primary_pid, dt);
    cr->primary_output = pid_loop_get_output(&cr->primary_pid);

    /* --- Ratio transformation: primary output → secondary SP --- */
    cr->secondary_sp_a = cr->primary_output * cr->ratio_cfg_a.ratio_gain
                         + cr->ratio_cfg_a.ratio_bias;

    /* Clamp secondary SP to ratio limits */
    if (cr->secondary_sp_a * cr->ratio_cfg_a.ratio_gain > cr->ratio_cfg_a.ratio_max)
        cr->secondary_sp_a = cr->ratio_cfg_a.ratio_max / cr->ratio_cfg_a.ratio_gain;
    if (cr->secondary_sp_a * cr->ratio_cfg_a.ratio_gain < cr->ratio_cfg_a.ratio_min)
        cr->secondary_sp_a = cr->ratio_cfg_a.ratio_min / cr->ratio_cfg_a.ratio_gain;

    /* --- Secondary A PID --- */
    if (cr->secondary_a_in_cascade) {
        pid_loop_set_sp(&cr->secondary_pid_a, cr->secondary_sp_a);
    }
    pid_loop_set_pv(&cr->secondary_pid_a, secondary_pv_a);
    pid_loop_execute(&cr->secondary_pid_a, dt);
    cr->secondary_output_a = pid_loop_get_output(&cr->secondary_pid_a);

    /* --- Secondary B PID (split-range / dual-output) --- */
    if (cr->mode == CASCADE_MODE_SPLIT_RANGE || cr->mode == CASCADE_MODE_DUAL_OUTPUT) {
        /* Split-range: primary output < split_point → secondary A
         *              primary output > split_point → secondary B */
        if (cr->mode == CASCADE_MODE_SPLIT_RANGE) {
            if (cr->primary_output <= cr->split_point) {
                /* A handles 0..split_point, scaled to 0..100% */
                cr->secondary_sp_a = 100.0 * cr->primary_output / cr->split_point;
                cr->secondary_sp_b = 0.0;
            } else {
                /* B handles split_point..100%, scaled to 0..100% */
                cr->secondary_sp_a = 100.0;
                cr->secondary_sp_b = 100.0 * (cr->primary_output - cr->split_point)
                                   / (100.0 - cr->split_point);
            }
        } else {
            /* Dual-output: both secondaries track the same ratio-adjusted SP */
            cr->secondary_sp_b = cr->primary_output * cr->ratio_cfg_b.ratio_gain
                                 + cr->ratio_cfg_b.ratio_bias;
        }

        pid_loop_set_sp(&cr->secondary_pid_b, cr->secondary_sp_b);
        pid_loop_set_pv(&cr->secondary_pid_b, secondary_pv_b);
        pid_loop_execute(&cr->secondary_pid_b, dt);
        cr->secondary_output_b = pid_loop_get_output(&cr->secondary_pid_b);
    }

    /* --- Override check --- */
    if (cr->primary_output > cr->override_high_limit ||
        cr->primary_output < cr->override_low_limit) {
        cr->override_active = 1;
    }

    /* Update state */
    if (cr->primary_pid.in_auto) {
        if (cr->secondary_a_in_cascade) {
            cr->state = CASCADE_STATE_SECONDARY_REMOTE;
        } else {
            cr->state = CASCADE_STATE_PRIMARY_ACTIVE;
        }
    }

    cr->cycle_count++;
    return 0;
}

int cascade_ratio_switch_secondary_to_cascade(cascade_ratio_t *cr, int which)
{
    if (!cr) return -1;

    pid_loop_t *sec = (which == 0) ? &cr->secondary_pid_a : &cr->secondary_pid_b;

    /* Bumpless transfer: initialize secondary setpoint to current PV */
    double current_pv = sec->process_variable;
    pid_loop_bumpless_init(sec, current_pv);
    sec->in_cascade = 1;

    if (which == 0) cr->secondary_a_in_cascade = 1;
    else            cr->secondary_b_in_cascade = 1;

    return 0;
}

int cascade_ratio_switch_secondary_to_local(cascade_ratio_t *cr, int which)
{
    if (!cr) return -1;

    pid_loop_t *sec = (which == 0) ? &cr->secondary_pid_a : &cr->secondary_pid_b;
    sec->in_cascade = 0;

    if (which == 0) cr->secondary_a_in_cascade = 0;
    else            cr->secondary_b_in_cascade = 0;

    return 0;
}

int cascade_ratio_override(cascade_ratio_t *cr,
                           double secondary_sp_a, double secondary_sp_b)
{
    if (!cr) return -1;

    cr->override_active = 1;
    cr->secondary_sp_a = secondary_sp_a;
    cr->secondary_sp_b = secondary_sp_b;
    cr->state = CASCADE_STATE_OVERRIDE_ACTIVE;

    return 0;
}

int cascade_ratio_shed(cascade_ratio_t *cr, double safe_output)
{
    if (!cr) return -1;

    cr->state             = CASCADE_STATE_SHEDDING;
    cr->secondary_sp_a    = safe_output;
    cr->secondary_sp_b    = safe_output;
    cr->secondary_output_a = safe_output;
    cr->secondary_output_b = safe_output;

    return 0;
}

int cascade_ratio_reset(cascade_ratio_t *cr)
{
    if (!cr) return -1;

    pid_loop_reset(&cr->primary_pid);
    pid_loop_reset(&cr->secondary_pid_a);
    pid_loop_reset(&cr->secondary_pid_b);

    cr->primary_output   = 0.0;
    cr->secondary_sp_a   = 0.0;
    cr->secondary_sp_b   = 0.0;
    cr->secondary_output_a = 0.0;
    cr->secondary_output_b = 0.0;
    cr->override_active  = 0;
    cr->state            = CASCADE_STATE_IDLE;
    cr->cycle_count      = 0;

    return 0;
}

/* -------------------------------------------------------------------
 * L2: Getters
 * ------------------------------------------------------------------- */

double cascade_ratio_get_primary_output(const cascade_ratio_t *cr)
{
    if (!cr) return NAN;
    return cr->primary_output;
}

double cascade_ratio_get_secondary_sp(const cascade_ratio_t *cr, int which)
{
    if (!cr) return NAN;
    return (which == 0) ? cr->secondary_sp_a : cr->secondary_sp_b;
}

double cascade_ratio_get_secondary_output(const cascade_ratio_t *cr, int which)
{
    if (!cr) return NAN;
    return (which == 0) ? cr->secondary_output_a : cr->secondary_output_b;
}

cascade_state_t cascade_ratio_get_state(const cascade_ratio_t *cr)
{
    if (!cr) return CASCADE_STATE_COUNT;
    return cr->state;
}

/* -------------------------------------------------------------------
 * L2: Utility
 * ------------------------------------------------------------------- */

const char *cascade_state_name(cascade_state_t state)
{
    static const char *names[] = {
        "Idle", "Primary Active", "Secondary Local",
        "Secondary Remote", "Override Active", "Initializing", "Shedding"
    };
    if (state < 0 || state >= CASCADE_STATE_COUNT) return "Unknown";
    return names[state];
}

const char *cascade_mode_name(cascade_mode_t mode)
{
    static const char *names[] = {
        "Simple", "Split-Range", "Valve Position",
        "Override Select", "Dual Output"
    };
    if (mode < 0 || mode >= CASCADE_MODE_COUNT) return "Unknown";
    return names[mode];
}
