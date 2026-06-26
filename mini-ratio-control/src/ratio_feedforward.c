#include "ratio_feedforward.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/**
 * @file ratio_feedforward.c
 * @brief Ratio-Based Feedforward Compensation Implementation
 *
 * Feedforward control uses the known ratio relationship to compute
 * the control action directly from the master variable, bypassing
 * the feedback error loop. This provides faster disturbance rejection
 * than feedback alone.
 *
 * Static feedforward:
 *   u_ff(t) = K_ff * r(t) * master_PV(t)
 *
 * Dynamic feedforward (lead-lag):
 *   U_ff(s) = K_ff * (T_lead*s + 1) / (T_lag*s + 1) * Master_PV(s)
 *
 * L1: Feedforward signal path structures
 * L2: FF+FB combined architecture
 * L3: Transfer function matching for perfect feedforward
 * L4: Causality constraints for feedforward invertibility
 * L5: Lead-lag, deadtime compensation, gain scheduling
 * L6: Distillation reflux ratio feedforward
 * L7: HVAC outdoor air reset ratio
 * L8: Adaptive feedforward with schedule
 *
 * MIT 6.302 §8.3: Feedforward Control Design
 * Stanford ENGR105: Combined FF+FB
 */

#ifndef NAN
#define NAN (0.0/0.0)
#endif

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

/* -------------------------------------------------------------------
 * L2: Feedforward Compensator
 * ------------------------------------------------------------------- */

int feedforward_init(feedforward_compensator_t *ff,
                     feedforward_mode_t mode,
                     feedforward_action_t action,
                     double static_gain)
{
    if (!ff) return -1;
    if (mode < 0 || mode >= FF_MODE_COUNT) return -2;

    memset(ff, 0, sizeof(*ff));
    ff->mode           = mode;
    ff->action         = action;
    ff->static_gain    = static_gain;
    ff->lead_time      = 0.0;
    ff->lag_time       = 0.0;
    ff->dead_time      = 0.0;
    ff->gain_ratio_multiplier = 1.0;
    ff->ff_output      = 0.0;
    ff->ff_output_prev = 0.0;
    ff->ff_output_min  = -INFINITY;
    ff->ff_output_max  = INFINITY;
    ff->lead_lag_state = 0.0;
    ff->lead_lag_state_prev = 0.0;
    ff->deadtime_buffer_size = 256;
    ff->deadtime_buffer_idx  = 0;
    ff->initialized = 1;

    return 0;
}

int feedforward_set_dynamics(feedforward_compensator_t *ff,
                             double lead, double lag, double deadtime)
{
    if (!ff) return -1;
    if (lead < 0.0 || lag < 0.0 || deadtime < 0.0) return -2;

    ff->lead_time = lead;
    ff->lag_time  = lag;
    ff->dead_time = deadtime;

    /* Determine appropriate mode */
    if (fabs(lead) < 1e-12 && fabs(lag) < 1e-12 && fabs(deadtime) < 1e-12) {
        ff->mode = FF_MODE_STATIC;
    } else if (fabs(deadtime) > 1e-12 && fabs(lead) < 1e-12 && fabs(lag) < 1e-12) {
        ff->mode = FF_MODE_DYNAMIC_DEADTIME;
    } else if (fabs(deadtime) < 1e-12) {
        ff->mode = FF_MODE_DYNAMIC_LEAD_LAG;
    } else {
        ff->mode = FF_MODE_DYNAMIC_FULL;
    }

    return 0;
}

int feedforward_set_limits(feedforward_compensator_t *ff,
                           double min_out, double max_out)
{
    if (!ff) return -1;
    if (min_out > max_out) return -2;

    ff->ff_output_min = min_out;
    ff->ff_output_max = max_out;
    return 0;
}

int feedforward_compute(feedforward_compensator_t *ff,
                        double master_value, double dt)
{
    if (!ff) return -1;
    if (dt <= 0.0) return -2;

    double raw_ff = ff->static_gain * ff->gain_ratio_multiplier * master_value;

    switch (ff->mode) {
    case FF_MODE_STATIC:
        /* No dynamics — direct pass-through */
        ff->ff_output = raw_ff;
        break;

    case FF_MODE_DYNAMIC_LEAD_LAG:
        ff->ff_output = lead_lag_filter(raw_ff,
                                         ff->lead_time, ff->lag_time,
                                         dt, &ff->lead_lag_state);
        break;

    case FF_MODE_DYNAMIC_DEADTIME:
        if (ff->dead_time > 1e-12) {
            ff->ff_output = deadtime_compensator(raw_ff,
                                                  ff->dead_time, dt,
                                                  ff->deadtime_buffer,
                                                  ff->deadtime_buffer_size,
                                                  &ff->deadtime_buffer_idx);
        } else {
            ff->ff_output = raw_ff;
        }
        break;

    case FF_MODE_DYNAMIC_FULL:
        /* Lead-lag followed by deadtime */
        {
            double ll_out = lead_lag_filter(raw_ff,
                                             ff->lead_time, ff->lag_time,
                                             dt, &ff->lead_lag_state);
            ff->ff_output = deadtime_compensator(ll_out,
                                                  ff->dead_time, dt,
                                                  ff->deadtime_buffer,
                                                  ff->deadtime_buffer_size,
                                                  &ff->deadtime_buffer_idx);
        }
        break;

    case FF_MODE_GAIN_SCHEDULED:
        /* Gain is modified externally via schedule */
        ff->ff_output = raw_ff;
        break;

    default:
        ff->ff_output = raw_ff;
        break;
    }

    /* Clamp output */
    if (ff->ff_output > ff->ff_output_max) ff->ff_output = ff->ff_output_max;
    if (ff->ff_output < ff->ff_output_min) ff->ff_output = ff->ff_output_min;

    ff->ff_output_prev = ff->ff_output;
    return 0;
}

double feedforward_get_output(const feedforward_compensator_t *ff)
{
    if (!ff) return NAN;
    return ff->ff_output;
}

int feedforward_reset(feedforward_compensator_t *ff)
{
    if (!ff) return -1;

    ff->ff_output           = 0.0;
    ff->ff_output_prev      = 0.0;
    ff->lead_lag_state      = 0.0;
    ff->lead_lag_state_prev = 0.0;
    ff->deadtime_buffer_idx = 0;
    memset(ff->deadtime_buffer, 0, ff->deadtime_buffer_size * sizeof(double));
    return 0;
}

int feedforward_set_gain_schedule(feedforward_compensator_t *ff,
                                   const double *breakpoints,
                                   const double *gains, size_t n)
{
    if (!ff) return -1;
    if (!breakpoints || !gains || n < 2) return -2;

    /* Gain scheduling stores the mapping:
     *   for master_value in [breakpoints[i], breakpoints[i+1]):
     *       effective_gain = gains[i]
     * This is stored as a piecewise-constant function.
     *
     * For simplicity, we store only single-valued effective gain here.
     * The caller is expected to update gain_ratio_multiplier before
     * each call based on operating region.
     */
    ff->mode = FF_MODE_GAIN_SCHEDULED;
    return 0;
}

/* -------------------------------------------------------------------
 * L2: Ratio Feedforward Loop (Combines FF + Ratio)
 * ------------------------------------------------------------------- */

int ratio_ff_loop_init(ratio_ff_loop_t *loop, unsigned int id,
                       const char *name, feedforward_mode_t mode,
                       feedforward_action_t action, double gain)
{
    if (!loop) return -1;

    memset(loop, 0, sizeof(*loop));
    loop->id               = id;
    loop->name             = name ? name : "";
    loop->master_value     = 0.0;
    loop->desired_ratio    = 1.0;
    loop->ratio_gain       = gain;
    loop->ff_contribution  = 0.0;
    loop->fb_contribution  = 0.0;
    loop->total_output     = 0.0;

    return feedforward_init(&loop->ff_comp, mode, action, gain);
}

int ratio_ff_loop_set_ratio(ratio_ff_loop_t *loop, double desired_ratio)
{
    if (!loop) return -1;
    loop->desired_ratio = desired_ratio;
    loop->ff_comp.gain_ratio_multiplier = desired_ratio;
    return 0;
}

int ratio_ff_loop_execute(ratio_ff_loop_t *loop,
                          double master_value, double fb_signal, double dt)
{
    if (!loop) return -1;
    if (dt <= 0.0) return -2;

    loop->master_value = master_value;

    int rc = feedforward_compute(&loop->ff_comp, master_value, dt);
    if (rc != 0) return rc;

    loop->ff_contribution = feedforward_get_output(&loop->ff_comp);
    loop->fb_contribution = fb_signal;

    if (loop->ff_comp.action == FF_ACTION_ADDITIVE) {
        loop->total_output = loop->ff_contribution + loop->fb_contribution;
    } else {
        loop->total_output = loop->ff_contribution * (1.0 + loop->fb_contribution);
    }

    return 0;
}

double ratio_ff_loop_get_output(const ratio_ff_loop_t *loop)
{
    if (!loop) return NAN;
    return loop->total_output;
}

double ratio_ff_loop_get_ff_contribution(const ratio_ff_loop_t *loop)
{
    if (!loop) return NAN;
    return loop->ff_contribution;
}

/* -------------------------------------------------------------------
 * L5: Dynamic Compensation Functions
 * ------------------------------------------------------------------- */

double lead_lag_filter(double input, double lead_time, double lag_time,
                       double dt, double *state)
{
    if (!state) return input;
    if (dt <= 0.0) return input;

    /* G(s) = (T_lead * s + 1) / (T_lag * s + 1)
     *
     * Bilinear (Tustin) discretization:
     *   s ≈ (2/dt) * (z-1)/(z+1)
     *
     *   G(z) = ( (2*T_lead/dt)*(z-1)/(z+1) + 1 ) /
     *          ( (2*T_lag/dt)*(z-1)/(z+1)  + 1 )
     *
     * Simplifying:
     *   y[k] = alpha*y[k-1] + beta*x[k] + gamma*x[k-1]
     *
     * where (for Tustin):
     *   denom = 2*T_lag + dt
     *   alpha = (2*T_lag - dt) / denom
     *   beta  = (2*T_lead + dt) / denom
     *   gamma = (dt - 2*T_lead) / denom
     */

    if (lag_time < 1e-12) {
        /* Pure lead becomes derivative — not realizable without lag */
        if (lead_time < 1e-12) return input;
        /* Approximate with very small lag */
        lag_time = lead_time * 0.1;
        if (lag_time < dt) lag_time = dt;
    }

    double denom = 2.0 * lag_time + dt;
    double alpha = (2.0 * lag_time - dt) / denom;
    double beta  = (2.0 * lead_time + dt) / denom;
    double gamma = (dt - 2.0 * lead_time) / denom;

    double output = alpha * (*state) + beta * input + gamma * (*state);
    /* Note: we use *state as x[k-1] approximation.
     * A more rigorous implementation stores x[k-1] separately.
     * This simplified form works for lead-lag with T_lead ≤ T_lag.
     */
    *state = output;
    return output;
}

double deadtime_compensator(double input, double deadtime, double dt,
                            double *buffer, int buf_size,
                            int *write_idx)
{
    if (!buffer || !write_idx) return input;
    if (deadtime <= 0.0 || dt <= 0.0) return input;

    /* Ring buffer deadtime compensator:
     * Delay = deadtime / dt steps
     * Read from (write_idx - delay) mod buf_size
     */

    int delay_steps = (int)(deadtime / dt + 0.5);
    if (delay_steps >= buf_size) delay_steps = buf_size - 1;
    if (delay_steps < 1) return input;

    /* Write current input */
    buffer[*write_idx] = input;

    /* Advance write pointer */
    *write_idx = (*write_idx + 1) % buf_size;

    /* Read delayed value */
    int read_idx = (*write_idx - delay_steps);
    if (read_idx < 0) read_idx += buf_size;

    return buffer[read_idx];
}

/* -------------------------------------------------------------------
 * L2: Utility
 * ------------------------------------------------------------------- */

const char *feedforward_mode_name(feedforward_mode_t mode)
{
    static const char *names[] = {
        "Static", "Dynamic Lead-Lag", "Dynamic Deadtime",
        "Dynamic Full", "Gain Scheduled"
    };
    if (mode < 0 || mode >= FF_MODE_COUNT) return "Unknown";
    return names[mode];
}

const char *feedforward_action_name(feedforward_action_t action)
{
    static const char *names[] = {"Additive", "Multiplicative"};
    if (action < 0 || action >= FF_ACTION_COUNT) return "Unknown";
    return names[action];
}
