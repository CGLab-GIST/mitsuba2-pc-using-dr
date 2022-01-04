#include <mitsuba/core/bbox.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/core/string.h>
#include <mitsuba/render/sensor.h>

NAMESPACE_BEGIN(mitsuba)


template <typename Float, typename Spectrum>
class Projector final : public Emitter<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Emitter, m_flags, m_needs_sample_3)
    MTS_IMPORT_TYPES(Texture)

    Projector(const Properties &props) : Base(props) {

        m_pos = props.point3f("position");
        // [roll pitch yaw]
        m_rot = props.vector3f("rotation");

        m_intensity = Texture::D65(props.float_("scale", 1));

        m_irradiance = props.texture<Texture>("irradiance");
        ScalarVector2i size = m_irradiance->resolution();
        m_x_fov = parse_fov(props, size.x() / (float) size.y());

        m_camera_to_sample = perspective_projection<Float>(
            size, size, 0, m_x_fov, 1e-4f, 1e4f);


        m_flags = +EmitterFlags::DeltaPosition;

        // 105% : 55
        // 115% : 65
        m_shift_y = props.float_("shift_y", 0);
        m_blur_size = props.float_("blur_size", 0);

    }

    /// TODO: Completely untested, need to revist once the particle tracer is merged.
    // Not used for now
    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                          const Point2f & /* spatial_sample */,
                                          const Point2f & direction_sample,
                                          Mask active) const override {
        NotImplementedError("Projector::sample_ray()");
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);
        return std::make_pair(
            Ray3f(), Spectrum());
    }

    // it : hitpoint
    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f & /*sample*/,
                     Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);

        // 1. Transform the reference point into the local coordinate system
        Transform4f trafo = Transform4f::look_at_2(m_pos, m_rot);
        Point it_local = trafo.inverse().transform_affine(it.p);

        // 2. Map to UV coordinates
        Point2f uv = head<2>(m_camera_to_sample * it_local);
        uv[1] += m_shift_y;

        active &= all(uv >= Point2f(0, 0) && uv <= Point2f(1, 1 )) && it_local.z() > 0;

        // 3. Query texture
        SurfaceInteraction3f it_query = zero<SurfaceInteraction3f>();
        it_query.wavelengths = it.wavelengths;
        it_query.uv = uv;
        UnpolarizedSpectrum spec = m_irradiance->eval(it_query, active);

        // 4. Prepare DirectionSample record for caller (MIS, etc.)
        DirectionSample3f ds;
        ds.p = m_pos;
        ds.n = trafo * Vector3f(0, 0, 1);
        ds.uv = uv;
        ds.time = it.time;
        ds.pdf = 1.f;
        ds.delta = true;
        ds.object   = this;

        ds.d = ds.p - it.p;
        Float dist_squared = squared_norm(ds.d);
        ds.dist = sqrt(dist_squared);
        ds.d *= rcp(ds.dist);

        // Scale so that irradiance at z=1 is correct
        spec *= math::Pi<Float> * m_intensity->eval(it_query, active) * smooth_profile(uv) *
                sqr(rcp(it_local.z())) / -dot(ds.n, ds.d);

        return { ds, unpolarized<Spectrum>(spec & active) };
    }

    Float pdf_direction(const Interaction3f &, const DirectionSample3f &,
                        Mask) const override {
        return 0.f;
    }

    Spectrum eval(const SurfaceInteraction3f & /*si*/,
                  Mask /*active*/) const override {
        return 0.f;
    }

    Float smooth_profile(const Point2f &uv) const {
        Float res_x(0);
        Float x = uv.x();
        res_x = select(x >= m_blur_size && x <= Float(1) - m_blur_size, Float(1), res_x);
        res_x = select(x < m_blur_size && x > Float(0), x / m_blur_size, res_x);
        res_x = select(x > Float(1) - m_blur_size && x < Float(1),
                     (1 - x) / m_blur_size, res_x);

        Float y = uv.y();
        Float res_y(0);
//        res_y = select(y >= m_blur_size - m_shift_y && y <= Float(1) - m_blur_size - m_shift_y, Float(1), res_y);
//        res_y = select(y < m_blur_size - m_shift_y && y > (-m_shift_y), (y+m_shift_y) / m_blur_size, res_y);
//        res_y = select(y > Float(1) - m_blur_size - m_shift_y && y < Float(1) - m_shift_y,
//                     (1 - y - m_shift_y) / m_blur_size, res_y);

        res_y = select(y >= m_blur_size && y <= Float(1) - m_blur_size, Float(1), res_y);
        res_y = select(y < m_blur_size && y > Float(0), y / m_blur_size, res_y);
        res_y = select(y > Float(1) - m_blur_size && y < Float(1),
                       (1 - y) / m_blur_size, res_y);
        
        return res_x * res_y;
    }

    ScalarBoundingBox3f bbox() const override {
        /* This emitter does not occupy any particular region
           of space, return an invalid bounding box */
        return ScalarBoundingBox3f();
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("irradiance", m_irradiance.get());
        callback->put_parameter("position", m_pos);
        callback->put_parameter("rotation", m_rot);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "Projector[" << std::endl
            << "  x_fov = " << m_x_fov << "," << std::endl
            << "  irradiance = " << string::indent(m_irradiance) << "," << std::endl
            << "  origin = " << m_pos << std::endl
            << "  dir = " << m_rot << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_irradiance;
    ref<Texture> m_intensity;
    Transform4f m_camera_to_sample;
    ScalarFloat m_x_fov;

    Point3f m_pos;
    Vector3f m_rot;

    Float m_blur_size;
    Float m_shift_y;
};


MTS_IMPLEMENT_CLASS_VARIANT(Projector, Emitter)
MTS_EXPORT_PLUGIN(Projector, "Projector emitter")
NAMESPACE_END(mitsuba)
