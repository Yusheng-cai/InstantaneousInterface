#include "Field.h"
#include "Mesh.h"
#include "tools/CommonTypes.h"


using Real =  CommonTypes::Real;

class Surface_nets{
	public:
		using Real = CommonTypes::Real;
		using Real3= CommonTypes::Real3;
		using INT3 = CommonTypes::index3;

		Surface_nets(){};

		void triangulate(const Field& field, Mesh& mesh, Real isoval = 0);

		int get_idx_from_ijk(const INT3& ijk);

    	// mapping from 3d coordinates to 1d
		INT3 get_ijk_from_idx(int idx);

	private:
		Real3 spacing_;
		INT3  N_;
};