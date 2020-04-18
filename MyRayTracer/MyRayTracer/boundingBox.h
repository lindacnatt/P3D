
#include "vector.h"
#include "ray.h"
#include "scene.h"

class AABB
{
public:
	Vector min, max;

	AABB(void);
	virtual ~AABB();
	AABB(const Vector& v0, const Vector& v1);
	AABB(const AABB& bbox);
	AABB operator= (const AABB& rhs);
	
	bool intercepts(const Ray& r, float& t, ShadeRec& sr);
	bool isInside(const Vector& p);
};