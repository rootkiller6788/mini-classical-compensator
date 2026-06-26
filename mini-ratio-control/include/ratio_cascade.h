#ifndef RATIO_CASCADE_H
#define RATIO_CASCADE_H

#include "ratio_control.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file ratio_cascade.h
 * @brief Cascade Ratio Control — multi-level ratio architectures
 *
 * Cascade control places one controller inside another, with the
 * inner (secondary) loop handling fast disturbances and the outer
 * (primary) loop handling slow setpoint tracking. Ratio cascade
 * adds ratio constraints between cascade levels.
 *
 * Architecture:
 *   Primary (master)  → produces demand signal
 *   Ratio station     → converts demand to secondary SP via ratio
 *   Secondary (slave) → tracks ratio-adjusted SP, rejects fast disturbances
 *
 * Common forms:
 *   1. Simple cascade: master_output → ratio → slave_SP
 *   2. Split-range cascade: master_output → ratio → two slave_SPs
 *   3. Valve position control: master_output → ratio → valve position SP
 *   4. Override cascade: high/low select between multiple ratio paths
 *
 * L1: Cascade loop definitions, hierarchy
 * L2: Cascade stability, inner loop faster than outer (N > 5 rule)
 * L3: Cascade transfer function decomposition
 * L4: Cascade integrity — outer loop bumpless transfer on inner mode change
 * L5: Anti-windup in cascade ratio loops
 * L6: Distillation column temperature/flow cascade
 * L7: Steam header pressure/flow cascade ratio in power plants
 * L8: Model-based cascade ratio with decoupling
 *
 * MIT 6.302 §9: Cascade Control
 * Cambridge 3F2: Cascade and Ratio
 * Purdue ME 675: Multivariable cascade
 */

typedef enum {
    CASCADE_STATE_IDLE = 0,
    CASCADE_STATE_PRIMARY_ACTIVE,
    CASCADE_STATE_SECONDARY_LOCAL,
    CASCADE_STATE_SECONDARY_REMOTE,
    CASCADE_STATE_OVERRIDE_ACTIVE,
    CASCADE_STATE_INITIALIZING,
    CASCADE_STATE_SHEDDING,
    CASCADE_STATE_COUNT
} cascade_state_t;

typedef enum {
    CASCADE_MODE_SIMPLE = 0,
    CASCADE_MODE_SPLIT_RANGE,
    CASCADE_MODE_VALVE_POSITION,
    CASCADE_MODE_OVERRIDE_SELECT,
    CASCADE_MODE_DUAL_OUTPUT,
    CASCADE_MODE_COUNT
} cascade_mode_t;

typedef struct {
    const char     *tag_name;
    double          setpoint;
    double          process_variable;
    double          output;
    double          output_min;
    double          output_max;
    double          kp, ki, kd;
    double          integral;
    double          error_prev;
    double          dt_accum;
    int             in_cascade;
    int             in_auto;
    double          antiwindup_limit;
    int             bumpless_active;
} pid_loop_t;

typedef struct {
    cascade_mode_t      mode;
    cascade_state_t     state;
    ratio_loop_t       *primary_ratio;
    pid_loop_t          primary_pid;
    pid_loop_t          secondary_pid_a;
    pid_loop_t          secondary_pid_b;
    double              primary_sp;
    double              primary_pv;
    double              primary_output;
    double              secondary_sp_a;
    double              secondary_sp_b;
    double              secondary_output_a;
    double              secondary_output_b;
    double              split_point;
    ratio_config_t      ratio_cfg_a;
    ratio_config_t      ratio_cfg_b;
    int                 secondary_a_in_cascade;
    int                 secondary_b_in_cascade;
    double              override_high_limit;
    double              override_low_limit;
    int                 override_active;
    double              init_target;
    unsigned long       cycle_count;
} cascade_ratio_t;

int  cascade_ratio_init(cascade_ratio_t *cr,
                        cascade_mode_t mode,
                        unsigned int primary_id,
                        const char *tag);
int  cascade_ratio_configure(cascade_ratio_t *cr,
                             const ratio_config_t *cfg_a,
                             const ratio_config_t *cfg_b);
int  cascade_ratio_tune_primary(cascade_ratio_t *cr,
                                double kp, double ki, double kd,
                                double out_min, double out_max);
int  cascade_ratio_tune_secondary(cascade_ratio_t *cr,
                                  int which,
                                  double kp, double ki, double kd,
                                  double out_min, double out_max);
int  cascade_ratio_set_primary_sp(cascade_ratio_t *cr, double sp);
int  cascade_ratio_execute(cascade_ratio_t *cr,
                           double primary_pv,
                           double secondary_pv_a,
                           double secondary_pv_b,
                           double dt);
int  cascade_ratio_switch_secondary_to_cascade(cascade_ratio_t *cr, int which);
int  cascade_ratio_switch_secondary_to_local(cascade_ratio_t *cr, int which);
int  cascade_ratio_override(cascade_ratio_t *cr,
                            double secondary_sp_a, double secondary_sp_b);
int  cascade_ratio_shed(cascade_ratio_t *cr, double safe_output);
int  cascade_ratio_reset(cascade_ratio_t *cr);

double cascade_ratio_get_primary_output(const cascade_ratio_t *cr);
double cascade_ratio_get_secondary_sp(const cascade_ratio_t *cr, int which);
double cascade_ratio_get_secondary_output(const cascade_ratio_t *cr, int which);
cascade_state_t cascade_ratio_get_state(const cascade_ratio_t *cr);

const char *cascade_state_name(cascade_state_t state);
const char *cascade_mode_name(cascade_mode_t mode);

int   pid_loop_init(pid_loop_t *pid, const char *tag,
                    double kp, double ki, double kd,
                    double out_min, double out_max);
int   pid_loop_set_sp(pid_loop_t *pid, double sp);
int   pid_loop_set_pv(pid_loop_t *pid, double pv);
int   pid_loop_execute(pid_loop_t *pid, double dt);
double pid_loop_get_output(const pid_loop_t *pid);
int   pid_loop_set_auto(pid_loop_t *pid, int auto_mode);
int   pid_loop_bumpless_init(pid_loop_t *pid, double current_output);
int   pid_loop_reset(pid_loop_t *pid);

#ifdef __cplusplus
}
#endif

#endif /* RATIO_CASCADE_H */
