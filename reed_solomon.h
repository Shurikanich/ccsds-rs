#ifndef INCLUDED_REED_SOLOMON_H
#define INCLUDED_REED_SOLOMON_H

#include <stdint.h>
#include <stdbool.h>



class reed_solomon
{
    private:
    public:
        reed_solomon();
        ~reed_solomon();

        void encode(uint8_t *data, bool use_dual_basis);
        int16_t decode(uint8_t *data, bool use_dual_basis);
};



#endif /* INCLUDED_REED_SOLOMON_H */
