#include "DensityField.h"

#include <string>

class AverageField:public DensityField
{
    public:
        AverageField(const DensityFieldInput& input);

        virtual void calculate() override;
        virtual void update() {};
    private:
};