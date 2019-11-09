#include "BSpline.h"
#include <cassert>
BSpline::BSpline(const std::vector<glm::vec3>& cpset) :
	m_control_points(cpset)
{
	recalc_knots();
}

BSpline::BSpline()
{
}

BSpline::~BSpline()
{
}

BSpline::BSpline(std::vector<glm::vec3>&& cpset) :
	m_control_points(std::move(cpset))
{
	recalc_knots();
}

void BSpline::appendControlPoint(const glm::vec3 & cp)
{
	m_control_points.push_back(cp);
	recalc_knots();
}

glm::vec3 BSpline::at(float t)
{
	if (m_control_points.size() - 1 < 3)
		return glm::vec3(0.0f);
	// clamp to [0, 1]
	t = glm::max(0.0f, glm::min(1.0f, t));
	// index of "active" knot interval
	int k = 0;
	for (int i = 3; i < m_knots.size() - 4; ++i)
	{
		k = i;
		if (t >= m_knots[i] && t < m_knots[i + 1])
		{			
			break;
		}
	}
	glm::vec3 d[4];
	for (int j = 0; j < 4; ++j)
		d[j] = m_control_points[j + k - 3];

	for (int r = 1; r <= 3; ++r)
	{
		for (int i = k - 3 + r; i <= k; ++i)
		{
			int j = i - (k - 3 + r);
			float a = (t - m_knots[i]) / (m_knots[i + 1 + 3 - r] - m_knots[i]);
			d[j] = (1.0f - a) * d[j] + a * d[j + 1];
		}
	}
	return d[0];
}

void BSpline::clear()
{
	m_control_points.clear();
}

void BSpline::recalc_knots()
{
	size_t m = m_control_points.size() - 1;
	constexpr size_t n = 3;
	if (n > m)
		return;
	float kvdelta = 1.0f / static_cast<float>(m - n + 1);
	m_knots.clear();
	// start padding
	for (size_t i = 0; i <= n; ++i)
		m_knots.push_back(0.0f);
	// >n, <= m
	for (size_t i = n + 1; i <= m; ++i)
		m_knots.push_back(static_cast<float>(i - n) * kvdelta);
	// end padding
	for (size_t i = m + 1; i <= m + n + 1; ++i)
		m_knots.push_back(1.0f);
}
