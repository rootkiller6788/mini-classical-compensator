#ifndef RATIO_FEEDFORWARD_H
#define RATIO_FEEDFORWARD_H

#include "ratio_control.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file ratio_feedforward.h
 * @brief Ratio-Based Feedforward Compensation
 *
 * Feedforward uses the ratio relationship to anticipate required
 * slave changes before the feedback controller sees error.
 *
 * Combined feedforward + feedback (FF+FB):
 *   u(t) = u_fb(t) + K_ff * r(t)
 * where u_fb = PID feedback output, u_ff = feedforward contribution.
 *
 * Static ratio FF:  u_ff = K_gain * master_PV * ratio
 * Dynamic ratio FF: u_ff = G_ff(s) * master_PV, where G_ff(s)
 *                    compensates for process dynamics differences.
 *
 * L1: Feedforward signal path definitions
 * L2: FF+FB architecture, static vs dynamic FF
 * L3: Transfer function matching for dynamic FF
 * L4: Perfect feedforward = process inverse (causality constraints)
 * L5: Lead-lag dynamic compensation, deadtime handling
 * L6: Distillation column feedforward ratio
 * L7: Ethylene cracker feed/steam ratio, HVAC outdoor reset
 * L8: Adaptive feedforward gain scheduling
 *
 * MIT 6.302 §8.3: Feedforward design
 * Stanford ENGR105: Feedforward + Feedback
 */

typedef enum {
    FF_MODE_STATIC = 0,
    FF_MODE_DYNAMIC_LEAD_LAG,
    FF_MODE_DYNAMIC_DEADTIME,
    FF_MODE_DYNAMIC_FULL,
    FF_MODE_GAIN_SCHEDULED,
    FF_MODE_COUNT
} feedforward_mode_t;

typedef enum {
    FF_ACTION_ADDITIVE = 0,
    FF_ACTION_MULTIPLICATIVE,
    FF_ACTION_COUNT
} feedforward_action_t;

typedef struct {
    feedforward_mode_t  mode;
    feedforward_action_t action;
    double              static_gain;
    double              lead_time;
    double              lag_time;
    double              dead_time;
    double              gain_ratio_multiplier;
    double              ff_output;
    double              ff_output_prev;
    double              ff_output_min;
    double              ff_output_max;
    double              deadtime_buffer[256];
    int                 deadtime_buffer_size;
    int                 deadtime_buffer_idx;
    double              lead_lag_state;
    double              lead_lag_state_prev;
    int                 initialized;
} feedforward_compensator_t;

typedef struct {
    unsigned int     id;
    const char      *name;
    double           master_value;
    double           desired_ratio;
    double           ratio_gain;
    double           ff_contribution;
    double           fb_contribution;
    double           total_output;
    feedforward_compensator_t ff_comp;
} ratio_ff_loop_t;

int  feedforward_init(feedforward_compensator_t *ff,
                      feedforward_mode_t mode,
                      feedforward_action_t action,
                      double static_gain);
int  feedforward_set_dynamics(feedforward_compensator_t *ff,
                              double lead, double lag, double deadtime);
int  feedforward_set_limits(feedforward_compensator_t *ff,
                            double min_out, double max_out);
int  feedforward_compute(feedforward_compensator_t *ff,
                         double master_value, double dt);
double feedforward_get_output(const feedforward_compensator_t *ff);
int  feedforward_reset(feedforward_compensator_t *ff);
int  feedforward_set_gain_schedule(feedforward_compensator_t *ff,
                                   const double *breakpoints,
                                   const double *gains, size_t n);

int  ratio_ff_loop_init(ratio_ff_loop_t *loop, unsigned int id,
                        const char *name, feedforward_mode_t mode,
                        feedforward_action_t action, double gain);
int  ratio_ff_loop_set_ratio(ratio_ff_loop_t *loop, double desired_ratio);
int  ratio_ff_loop_execute(ratio_ff_loop_t *loop,
                           double master_value, double fb_signal, double dt);
double ratio_ff_loop_get_output(const ratio_ff_loop_t *loop);
double ratio_ff_loop_get_ff_contribution(const ratio_ff_loop_t *loop);

const char *feedforward_mode_name(feedforward_mode_t mode);
const char *feedforward_action_name(feedforward_action_t action);

double lead_lag_filter(double input, double lead_time, double lag_time,
                       double dt, double *state);
double deadtime_compensator(double input, double deadtime, double dt,
                            double *buffer, int buf_size,
                            int *write_idx);

#ifdef __cplusplus
}
#endif

#endif /* RATIO_FEEDFORWARD_H */
