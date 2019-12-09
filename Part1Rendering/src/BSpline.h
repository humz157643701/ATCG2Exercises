#ifndef _B_SPLINE_H_
#define _B_SPLINE_H_

#include <libheaders.h>
#include <vector>

// simple 3d cubic b-spline implementation with uniformly spaced knots
class BSpline
{
public:
	BSpline();
	~BSpline();

	BSpline(const BSpline& other) = default;
	BSpline(BSpline&& other) = default;

	BSpline& operator=(const BSpline& other) = default;
	BSpline& operator=(BSpline&& other) = default;

	BSpline(const std::vector<glm::vec3>& cpset);
	BSpline(std::vector<glm::vec3>&& cpset);
	void appendControlPoint(const glm::vec3& cp);
	glm::vec3 at(float t);
	void clear();

private:
	void recalc_knots();

	std::vector<glm::vec3> m_control_points;
	std::vector<float> m_knots;
};







#endif