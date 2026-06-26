#ifndef RATIO_CONTROL_H
#define RATIO_CONTROL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    RATIO_MODE_MASTER_SLAVE  = 0,
    RATIO_MODE_CROSS_LIMITING,
    RATIO_MODE_BLENDING,
    RATIO_MODE_FEEDFORWARD,
    RATIO_MODE_CASCADE_RATIO,
    RATIO_MODE_ADAPTIVE,
    RATIO_MODE_SPLIT_RANGE,
    RATIO_MODE_COUNT
} ratio_mode_t;

typedef enum {
    RATIO_STATION_INTERNAL = 0,
    RATIO_STATION_EXTERNAL,
    RATIO_STATION_COMPUTED
} ratio_station_type_t;

typedef enum {
    RATIO_ACTION_MULTIPLY = 0,
    RATIO_ACTION_DIVIDE,
    RATIO_ACTION_PERCENT,
    RATIO_ACTION_CHARACTERIZED
} ratio_action_t;

typedef enum {
    RATIO_INTEGRITY_OK = 0,
    RATIO_INTEGRITY_MASTER_BAD,
    RATIO_INTEGRITY_SLAVE_BAD,
    RATIO_INTEGRITY_RATIO_CLAMPED,
    RATIO_INTEGRITY_SLAVE_SATURATED,
    RATIO_INTEGRITY_CROSS_LIMITED,
    RATIO_INTEGRITY_MASS_BAL_VIOLATED,
    RATIO_INTEGRITY_SLAVE_TRACKING_ERR,
    RATIO_INTEGRITY_COUNT
} ratio_integrity_t;

typedef enum {
    RATIO_DIRECT_ACTING = 0,
    RATIO_REVERSE_ACTING
} ratio_direction_t;

typedef struct {
    ratio_mode_t        mode;
    ratio_station_type_t station_type;
    ratio_action_t      action;
    ratio_direction_t   direction;
    double              ratio_gain;
    double              ratio_bias;
    double              ratio_min;
    double              ratio_max;
    double              rate_limit_up;
    double              rate_limit_down;
    double              master_span;
    double              slave_span;
    double              master_zero;
    double              slave_zero;
    int                 filter_enabled;
    double              filter_tau;
    double              deadband;
    double              slave_sp_min;
    double              slave_sp_max;
    int                 shed_on_master_fail;
    double              shed_slave_sp;
} ratio_config_t;

typedef struct {
    double              master_pv;
    double              master_pv_filt;
    double              slave_pv;
    double              slave_pv_filt;
    double              slave_sp;
    double              slave_sp_prev;
    double              actual_ratio;
    double              ratio_deviation;
    double              slave_error;
    double              master_rate;
    double              master_prev;
    ratio_integrity_t   integrity;
    int                 master_valid;
    int                 slave_valid;
    int                 shed_active;
    unsigned long       cycle_count;
    double              total_master_flow;
    double              total_slave_flow;
    double              integrator_dt;
} ratio_state_t;

typedef struct {
    double              input;
    double              ratio_value;
} ratio_char_point_t;

typedef struct {
    ratio_char_point_t *points;
    size_t              num_points;
    double              extrap_low;
    double              extrap_high;
} ratio_char_table_t;

typedef struct {
    unsigned int        loop_id;
    const char         *tag;
    const char         *description;
    const char         *master_tag;
    const char         *slave_tag;
    const char         *eng_units_master;
    const char         *eng_units_slave;
    ratio_config_t      config;
    ratio_state_t       state;
    ratio_char_table_t *char_table;
    void               *user_data;
} ratio_loop_t;

int  ratio_init(ratio_loop_t *loop, unsigned int loop_id,
                const char *tag, const char *master_tag, const char *slave_tag);
int  ratio_configure(ratio_loop_t *loop, ratio_mode_t mode,
                     ratio_action_t action, double gain, double bias,
                     double rmin, double rmax);
int  ratio_set_master(ratio_loop_t *loop, double master_pv, double dt);
int  ratio_set_slave(ratio_loop_t *loop, double slave_pv);
double ratio_get_slave_sp(const ratio_loop_t *loop);
double ratio_get_actual(const ratio_loop_t *loop);
ratio_integrity_t ratio_get_integrity(const ratio_loop_t *loop);
int  ratio_set_gain(ratio_loop_t *loop, double new_gain, double dt);
int  ratio_shed(ratio_loop_t *loop, double safe_sp);
int  ratio_reset_state(ratio_loop_t *loop);

int  ratio_char_table_init(ratio_char_table_t *table,
                           const ratio_char_point_t *points, size_t n);
double ratio_char_lookup(const ratio_char_table_t *table, double input);
void ratio_char_table_free(ratio_char_table_t *table);

const char *ratio_mode_name(ratio_mode_t mode);
const char *ratio_integrity_name(ratio_integrity_t status);
const char *ratio_action_name(ratio_action_t action);
int         ratio_validate_config(const ratio_config_t *cfg);
double      ratio_check_mass_balance(const double *fractions, size_t n, double tol);
double      ratio_optimal_blend(double cost_a, double cost_b,
                                double rmin, double rmax);

#ifdef __cplusplus
}
#endif

#endif /* RATIO_CONTROL_H */
