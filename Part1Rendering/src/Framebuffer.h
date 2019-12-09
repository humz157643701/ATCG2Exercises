#ifndef _FRAMEBUFFER_H_
#define _FRAMEBUFFER_H_

class RenderBuffer
{
public:
	RenderBuffer() :
		rbid(0)
	{}
	RenderBuffer(GLuint id) :
		rbid(id)
	{}
	GLuint rbid;

	~RenderBuffer()
	{
		if (rbid)
			glDeleteRenderbuffers(1, &rbid);
	}
};

enum class RenderTargetType
{
	//TODO: add support for commented target types!
	Empty,
	RenderBuffer,
	//RenderBufferMS,
	Texture2D,
	//Texture2DMS,
	TextureCube
};

class RenderTarget
{
public:
	RenderTarget() :
		attachment(GL_INVALID_ENUM),
		rb(),
		type(RenderTargetType::Empty),
		width(0),
		height(0),
		miplevels(0),
		active(false)
	{}
	~RenderTarget()
	{}
	RenderTarget(std::shared_ptr<Texture2D> _tex, GLenum _attachment, GLsizei _width, GLsizei _height, GLsizei _miplevels = 1) :
		tex(_tex),
		attachment(_attachment),
		type(RenderTargetType::Texture2D),
		width(_width),
		height(_height),
		miplevels(_miplevels),
		active(true)
	{}
	RenderTarget(std::shared_ptr<RenderBuffer> _rb, GLenum _attachment, GLsizei _width, GLsizei _height) :
		rb(_rb),
		attachment(_attachment),
		type(RenderTargetType::RenderBuffer),
		width(_width),
		height(_height),
		miplevels(1),
		active(true)
	{}
	RenderTarget(std::shared_ptr<TextureCube> _ctex, GLenum _attachment, GLsizei _width, GLsizei _height, GLsizei _miplevels = 1) :
		ctex(_ctex),
		attachment(_attachment),
		type(RenderTargetType::TextureCube),
		width(_width),
		height(_height),
		miplevels(_miplevels),
		active(true)
	{}
	RenderTarget(const RenderTarget& other) 
	{
		if (other.type == RenderTargetType::Empty)
		{
			type = RenderTargetType::Empty;
			attachment = GL_INVALID_ENUM;
			width = other.width;
			height = other.height;
			miplevels = other.miplevels;
			active = other.active;
		}
		if (other.type == RenderTargetType::Texture2D)
		{
			attachment = other.attachment;
			tex = other.tex;
			type = RenderTargetType::Texture2D;
			width = other.width;
			height = other.height;
			miplevels = other.miplevels;
			active = other.active;
		}
		else if (other.type == RenderTargetType::RenderBuffer)
		{
			attachment = other.attachment;
			rb = other.rb;
			type = RenderTargetType::RenderBuffer;
			width = other.width;
			height = other.height;
			miplevels = other.miplevels;
			active = other.active;
		}
		else if (other.type == RenderTargetType::TextureCube)
		{
			attachment = other.attachment;
			ctex = other.ctex;
			type = RenderTargetType::TextureCube;
			width = other.width;
			height = other.height;
			miplevels = other.miplevels;
			active = other.active;
		}
		else //cube map follows
		{
			
		}
	}
	RenderTarget& operator=(const RenderTarget& other)
	{
		if (this == &other)
			return *this;
		if (other.type == RenderTargetType::Empty)
		{
			type = RenderTargetType::Empty;
			attachment = GL_INVALID_ENUM;
			width = other.width;
			height = other.height;
			miplevels = other.miplevels;
			active = other.active;
		}
		if (other.type == RenderTargetType::Texture2D)
		{
			attachment = other.attachment;
			tex = other.tex;
			type = RenderTargetType::Texture2D;
			width = other.width;
			height = other.height;
			miplevels = other.miplevels;
			active = other.active;
		}
		else if (other.type == RenderTargetType::RenderBuffer)
		{
			attachment = other.attachment;
			rb = other.rb;
			type = RenderTargetType::RenderBuffer;
			width = other.width;
			height = other.height;
			miplevels = other.miplevels;
			active = other.active;
		}
		else if (other.type == RenderTargetType::TextureCube)
		{
			attachment = other.attachment;
			ctex = other.ctex;
			type = RenderTargetType::TextureCube;
			width = other.width;
			height = other.height;
			miplevels = other.miplevels;
			active = other.active;
		}
		else //cube map follows
		{

		}
		return *this;
	}
	RenderTargetType type;
	std::shared_ptr<Texture> tex;
	GLsizei width;
	GLsizei height;
	GLsizei miplevels;
	GLenum attachment;
	bool active;
};

struct RenderTargetSet
{
	std::vector<RenderTarget> colorTargets;
	RenderTarget depthTarget;
};

class FrameBuffer
{
public:
	FrameBuffer() : 
		fbo(0),
		state(GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT)
	{}
	FrameBuffer(GLuint _fbo) :
		fbo(_fbo),
		state(GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT)
	{
		
	}
	~FrameBuffer()
	{
		if (fbo)
			glDeleteFramebuffers(1, &fbo);
	}
	void bind(GLenum target);
	void unbind(GLenum target);
	bool selectColorTargetMipmapLevel(size_t colorTargetIndex , GLint level);
	bool selectDepthTargetMipmapLevel(GLint level);
	bool attachRenderTargetSet(const RenderTargetSet& rtset);
	GLuint fbo;
	/*std::vector<RenderTarget> colorTargets;
	RenderTarget depthTarget;*/
	RenderTargetSet rtset;
	GLenum state;
	bool isComplete()
	{
		return state == GL_FRAMEBUFFER_COMPLETE;
	}
};

struct RenderTargetDesc
{
	RenderTargetDesc() :
		width(0),
		height(0),
		internalformat(GL_INVALID_ENUM),
		attachment(GL_INVALID_ENUM),
		type(RenderTargetType::Empty),
		samples(0),
		miplevels(0)
	{}

	RenderTargetDesc(GLsizei _width,
					 GLsizei _height,
					 GLenum _internalformat,
					 GLenum _attachment,
					 RenderTargetType _type,
					 int _miplevels = 1,
					 int samples = 0) :
		width(_width),
		height(_height),
		internalformat(_internalformat),
		attachment(_attachment),
		type(_type),
		samples(samples),
		miplevels(_miplevels)
	{}

	~RenderTargetDesc()
	{}

	GLsizei width;
	GLsizei height;
	GLenum internalformat;
	GLenum attachment;
	RenderTargetType type;
	int samples;
	int miplevels;
};

struct FrameBufferDesc
{
	std::vector<RenderTargetDesc> colorTargets;
	RenderTargetDesc depthTarget;
};

#endif