#ifndef TTV_ERRORS_H_
#define TTV_ERRORS_H_

/* Structure representing the state of a function which may be caused to return
 * early.
 */
typedef enum {
    STATUS_OK = 0,           // Convention for 0 to mean ok
    STATUS_ARGUMENT_ERROR,   // ValueError('Input flag must be 0, 1, or 2.')
    STATUS_MEMORY_ERROR,     // MemoryError('Not enough memory allocated for Transit
                             //   structure: more events triggering as transits than
                             //   expected. Possibily indicative of larger problem.')
    STATUS_HYPERBOLIC_ORBIT, // ValueError('Hyperbolic orbit.')
    STATUS_TOO_MANY_PLANETS, // ValueError('Too many planets.')
    STATUS_NON_CONVERGING,   // ValueError('Kepler step not converging in MAX_ITER.
                             //   Likely need a smaller dt.')
} status_t;

#endif // TTV_ERRORS_H_
