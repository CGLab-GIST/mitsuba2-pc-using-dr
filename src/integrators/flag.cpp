#include <enoki/stl.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/core/properties.h>

NAMESPACE_BEGIN(mitsuba)


template <typename Float, typename Spectrum>
class FlagIntegrator : public SamplingIntegrator<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(SamplingIntegrator, m_hide_emitters)
    MTS_IMPORT_TYPES(Scene, Sampler, Medium, Emitter, EmitterPtr, BSDF, BSDFPtr)

    // =============================================================
    //! @{ \name Constructors
    // =============================================================
    FlagIntegrator(const Properties &props) : Base(props) {
        if (props.has_property("shading_samples")
            && (props.has_property("emitter_samples") ||
                props.has_property("bsdf_samples"))) {
            Throw("Cannot specify both 'shading_samples' and"
                  " ('emitter_samples' and/or 'bsdf_samples').");
        }

        /// Number of shading samples -- this parameter is a shorthand notation
        /// to set both 'emitter_samples' and 'bsdf_samples' at the same time
        size_t shading_samples = props.size_("shading_samples", 1);

        /// Number of samples to take using the emitter sampling technique
        m_emitter_samples = props.size_("emitter_samples", shading_samples);

        /// Number of samples to take using the BSDF sampling technique
        m_bsdf_samples = props.size_("bsdf_samples", shading_samples);

        if (m_emitter_samples + m_bsdf_samples == 0)
            Throw("Must have at least 1 BSDF or emitter sample!");

        size_t sum    = m_emitter_samples + m_bsdf_samples;
    }

    std::pair<Spectrum, Mask> sample(const Scene *scene,
                                     Sampler *sampler,
                                     const RayDifferential3f &ray,
                                     const Medium * /* medium */,
                                     Float * /* aovs */,
                                     Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);

        SurfaceInteraction3f si = scene->ray_intersect(ray, active);
        Mask valid_ray = true;

        Spectrum result(0.f);

        // ----------------------- Visible emitters -----------------------
        if (!m_hide_emitters) {
            EmitterPtr emitter_vis = si.emitter(scene, active);
            if (any_or<true>(neq(emitter_vis, nullptr)))
                result += emitter_vis->eval(si, active);
        }

        active &= si.is_valid();
        if (none_or<false>(active))
            return { result, valid_ray };

        // ----------------------- Emitter sampling -----------------------

        BSDFContext ctx;
        BSDFPtr bsdf = si.bsdf(ray);
        Mask sample_emitter = active && has_flag(bsdf->flags(), BSDFFlags::Smooth);

        if (any_or<true>(sample_emitter)) {
            for (size_t i = 0; i < m_emitter_samples; ++i) {
                Mask active_e = sample_emitter;
                DirectionSample3f ds;
                Spectrum emitter_val;
                std::tie(ds, emitter_val) = scene->sample_emitter_direction(
                    si, sampler->next_2d(active_e), true, active_e);
                active_e &= neq(ds.pdf, 0.f);
                if (none_or<false>(active_e))
                    continue;

                emitter_val[emitter_val>1] = 1;
                
                result[active_e] += emitter_val;
            }
        }

        return { result, Mask(true) };
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "FlagIntegrator[" << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    size_t m_emitter_samples;
    size_t m_bsdf_samples;
};

MTS_IMPLEMENT_CLASS_VARIANT(FlagIntegrator, SamplingIntegrator)
MTS_EXPORT_PLUGIN(FlagIntegrator, "Direct integrator");
NAMESPACE_END(mitsuba)
