set_project("drift-retained-face-candidate-checks")
set_languages("c++17")
add_rules("mode.release")
add_requires("eigen", "tbb")
add_requires("suitesparse", {configs = {graphblas = false}})
add_defines("NOMINMAX", "_USE_MATH_DEFINES")
add_cxflags("/utf-8", "/bigobj", {tools = "cl"})
for _, name in ipairs({"curved", "preparation", "observer"}) do
    target("test_retained_" .. name)
        set_kind("binary")
        add_files("../../tests/test_retained_" .. name .. ".cpp")
        add_includedirs("../../include", "../../../core/include")
        add_packages("eigen", "tbb", "suitesparse")
        add_defines("VF_HAVE_CHOLMOD")
        set_targetdir(path.join(os.scriptdir(), "../../../results/local-runs/retained_checks_20261113/bin"))
        on_load(function (target)
            local pkg = target:pkg("suitesparse")
            if pkg and pkg:installdir() then
                target:add("includedirs", path.join(pkg:installdir(), "include", "suitesparse"))
            end
        end)
    target_end()
end
