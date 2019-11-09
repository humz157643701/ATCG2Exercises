#include "Camera.h"

Camera::Camera() :
	m_width(800),
	m_height(600),
	m_fovy(90.0f),
	m_near(0.1f),
	m_far(100.0f),
	m_aspectratio(800.0f / 600.0f),
	m_viewDirty(true),
	m_projDirty(true)
{}

Camera::Camera(int width, int height, float fovy, float near, float far) :
	m_width(static_cast<float>(width)),
	m_height(static_cast<float>(height)),
	m_fovy(fovy),
	m_near(near),
	m_far(far),
	m_aspectratio(static_cast<float>(width) / static_cast<float>(height)),
	m_viewDirty(true),
	m_projDirty(true)
{}

Camera::~Camera()
{}

const glm::mat4 & Camera::getViewMatrix()
{
	if (m_viewDirty)
	{
		auto wp = transform.getWorldPosition();
		m_viewMatrix = glm::lookAt(wp, wp - transform.getWorldZAxis(), transform.getWorldYAxis());
		m_viewDirty = false;
	}
	return m_viewMatrix;
}

const glm::mat4 & Camera::getProjectionMatrix()
{
	if (m_projDirty)
	{
		m_projectionMatrix = glm::perspective(m_fovy, m_aspectratio, m_near, m_far);
		m_projDirty = false;
	}
	return m_projectionMatrix;
}

bool Camera::bind(ShaderProgram* shader, const std::string & instanceName)
{
	bool success = true;
	success = shader->setUniform((instanceName + ".viewMatrix").c_str(), getViewMatrix(), false) && success;
	success = shader->setUniform((instanceName + ".projectionMatrix").c_str(), getProjectionMatrix(), false) && success;
	success = shader->setUniform((instanceName + ".worldPos").c_str(), transform.getWorldPosition()) && success;
	success = shader->setUniform((instanceName + ".renderWidth").c_str(), m_width) && success;
	success = shader->setUniform((instanceName + ".renderHeight").c_str(), m_height) && success;
	success = shader->setUniform((instanceName + ".nearPlane").c_str(), m_near) && success;
	success = shader->setUniform((instanceName + ".farPlane").c_str(), m_far) && success;
	return success;
}

void Camera::forward(float amount)
{
	transform.translateLocal(glm::vec3(0.0f, 0.0f, -amount));
	m_viewDirty = true;
}

void Camera::backward(float amount)
{
	transform.translateLocal(glm::vec3(0.0f, 0.0f, amount));
	m_viewDirty = true;
}

void Camera::left(float amount)
{
	transform.translateLocal(glm::vec3(-amount, 0.0f, 0.0f));
	m_viewDirty = true;
}

void Camera::right(float amount)
{
	transform.translateLocal(glm::vec3(amount, 0.0f, 0.0f));
	m_viewDirty = true;
}

void Camera::up(float amount)
{
	translate(glm::vec3(0.0f, amount, 0.0f));
	m_viewDirty = true;
}

void Camera::down(float amount)
{
	translate(glm::vec3(0.0f, -amount, 0.0f));
	m_viewDirty = true;
}

void Camera::rotateView(float pitch, float yaw)
{
	glm::quat pitchrot(glm::angleAxis(pitch, glm::vec3(1.0f, 0.0f, 0.0f)));
	glm::quat yawrot(glm::angleAxis(yaw, glm::vec3(0.0f, 1.0f, 0.0f)));
	transform.rotateLocal(pitchrot);
	transform.rotate(yawrot);
	m_viewDirty = true;
}

void Camera::setViewport(int width, int height)
{
	m_width = static_cast<float>(width);
	m_height = static_cast<float>(height);
	m_aspectratio = m_width / m_height;
	m_projDirty = true;
}

void Camera::setFrustumBounds(float near, float far)
{
	m_near = near;
	m_far = far;
	m_projDirty = true;
}

void Camera::setFieldOfView(float fovy)
{
	m_fovy = fovy;
	m_projDirty = true;
}

void Camera::setParameters(int width, int height, float fovy, float near, float far)
{
	setViewport(width, height);
	setFrustumBounds(near, far);
	setFieldOfView(fovy);
}

float Camera::getWidth()
{
	return m_width;
}

float Camera::getHeight()
{
	return m_height;
}

float Camera::getNear()
{
	return m_near;
}

float Camera::getFar()
{
	return m_far;
}

float Camera::getFieldOfView()
{
	return m_fovy;
}

float Camera::getAspectRatio()
{
	return m_aspectratio;
}

void Camera::setViewDirty()
{
	m_viewDirty = true;
}
