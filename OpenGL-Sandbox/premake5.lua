project "OpenGL-Sandbox"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++17"
	staticruntime "on"

	targetdir ("../bin/" .. outputdir .. "/%{prj.name}")
	objdir ("../bin-int/" .. outputdir .. "/%{prj.name}")

	files
	{
		"src/**.h",
		"src/**.cpp",
    "src/**.hpp"
	}

  libdirs 
  {
    "lib"
  }

	includedirs
	{
		"../PD-Core/src",
    "../OpenGL-Core/vendor/gltf",
		"../OpenGL-Core/vendor/spdlog/include",
		"../OpenGL-Core/src",
		"../OpenGL-Core/vendor",
		"../OpenGL-Core/%{IncludeDir.glm}",
		"../OpenGL-Core/%{IncludeDir.Glad}",
		"../OpenGL-Core/%{IncludeDir.ImGui}",
		"../OpenGL-Core/%{IncludeDir.assimp}/include",
		"%{IncludeDir.assimp}/contrib/irrXML",
		"%{IncludeDir.assimp}/contrib/zlib",
		"%{IncludeDir.assimp}/contrib/rapidjson/include",
    "../OpenGL-Core/%{IncludeDir.eigen}/Eigen"
	}

  defines 
  {
		-- Only PDSolver's own OpenMP loop (m_parallel_local_step) uses threads.
		"EIGEN_DONT_PARALLELIZE",
		-- "SWIG",
		"ASSIMP_BUILD_NO_OWN_ZLIB",

		"ASSIMP_BUILD_NO_X_IMPORTER",
		"ASSIMP_BUILD_NO_3DS_IMPORTER",
		"ASSIMP_BUILD_NO_MD3_IMPORTER",
		"ASSIMP_BUILD_NO_MDL_IMPORTER",
		"ASSIMP_BUILD_NO_MD2_IMPORTER",
		-- "ASSIMP_BUILD_NO_PLY_IMPORTER",
		"ASSIMP_BUILD_NO_ASE_IMPORTER",
		"ASSIMP_BUILD_NO_OBJ_IMPORTER",
		"ASSIMP_BUILD_NO_AMF_IMPORTER",
		"ASSIMP_BUILD_NO_HMP_IMPORTER",
		"ASSIMP_BUILD_NO_SMD_IMPORTER",
		"ASSIMP_BUILD_NO_MDC_IMPORTER",
		"ASSIMP_BUILD_NO_MD5_IMPORTER",
		"ASSIMP_BUILD_NO_STL_IMPORTER",
		"ASSIMP_BUILD_NO_LWO_IMPORTER",
		"ASSIMP_BUILD_NO_DXF_IMPORTER",
		"ASSIMP_BUILD_NO_NFF_IMPORTER",
		"ASSIMP_BUILD_NO_RAW_IMPORTER",
		"ASSIMP_BUILD_NO_OFF_IMPORTER",
		"ASSIMP_BUILD_NO_AC_IMPORTER",
		"ASSIMP_BUILD_NO_BVH_IMPORTER",
		"ASSIMP_BUILD_NO_IRRMESH_IMPORTER",
		"ASSIMP_BUILD_NO_IRR_IMPORTER",
		"ASSIMP_BUILD_NO_Q3D_IMPORTER",
		"ASSIMP_BUILD_NO_B3D_IMPORTER",
		-- "ASSIMP_BUILD_NO_COLLADA_IMPORTER",
		"ASSIMP_BUILD_NO_TERRAGEN_IMPORTER",
		"ASSIMP_BUILD_NO_CSM_IMPORTER",
		"ASSIMP_BUILD_NO_3D_IMPORTER",
		"ASSIMP_BUILD_NO_LWS_IMPORTER",
		"ASSIMP_BUILD_NO_OGRE_IMPORTER",
		"ASSIMP_BUILD_NO_OPENGEX_IMPORTER",
		"ASSIMP_BUILD_NO_MS3D_IMPORTER",
		"ASSIMP_BUILD_NO_COB_IMPORTER",
		"ASSIMP_BUILD_NO_BLEND_IMPORTER",
		"ASSIMP_BUILD_NO_Q3BSP_IMPORTER",
		"ASSIMP_BUILD_NO_NDO_IMPORTER",
		"ASSIMP_BUILD_NO_IFC_IMPORTER",
		"ASSIMP_BUILD_NO_XGL_IMPORTER",
		"ASSIMP_BUILD_NO_FBX_IMPORTER",
		"ASSIMP_BUILD_NO_ASSBIN_IMPORTER",
		-- "ASSIMP_BUILD_NO_GLTF_IMPORTER",
		"ASSIMP_BUILD_NO_C4D_IMPORTER",
		"ASSIMP_BUILD_NO_3MF_IMPORTER",
		"ASSIMP_BUILD_NO_X3D_IMPORTER",
		"ASSIMP_BUILD_NO_MMD_IMPORTER",
		
		"ASSIMP_BUILD_NO_STEP_EXPORTER",
		"ASSIMP_BUILD_NO_SIB_IMPORTER",

		-- "ASSIMP_BUILD_NO_MAKELEFTHANDED_PROCESS",
		-- "ASSIMP_BUILD_NO_FLIPUVS_PROCESS",
		-- "ASSIMP_BUILD_NO_FLIPWINDINGORDER_PROCESS",
		-- "ASSIMP_BUILD_NO_CALCTANGENTS_PROCESS",
		"ASSIMP_BUILD_NO_JOINVERTICES_PROCESS",
		-- "ASSIMP_BUILD_NO_TRIANGULATE_PROCESS",
		"ASSIMP_BUILD_NO_GENFACENORMALS_PROCESS",
		-- "ASSIMP_BUILD_NO_GENVERTEXNORMALS_PROCESS",
		"ASSIMP_BUILD_NO_REMOVEVC_PROCESS",
		"ASSIMP_BUILD_NO_SPLITLARGEMESHES_PROCESS",
		"ASSIMP_BUILD_NO_PRETRANSFORMVERTICES_PROCESS",
		"ASSIMP_BUILD_NO_LIMITBONEWEIGHTS_PROCESS",
		-- "ASSIMP_BUILD_NO_VALIDATEDS_PROCESS",
		"ASSIMP_BUILD_NO_IMPROVECACHELOCALITY_PROCESS",
		"ASSIMP_BUILD_NO_FIXINFACINGNORMALS_PROCESS",
		"ASSIMP_BUILD_NO_REMOVE_REDUNDANTMATERIALS_PROCESS",
		"ASSIMP_BUILD_NO_FINDINVALIDDATA_PROCESS",
		"ASSIMP_BUILD_NO_FINDDEGENERATES_PROCESS",
		"ASSIMP_BUILD_NO_SORTBYPTYPE_PROCESS",
		"ASSIMP_BUILD_NO_GENUVCOORDS_PROCESS",
		"ASSIMP_BUILD_NO_TRANSFORMTEXCOORDS_PROCESS",
		"ASSIMP_BUILD_NO_FINDINSTANCES_PROCESS",
		"ASSIMP_BUILD_NO_OPTIMIZEMESHES_PROCESS",
		"ASSIMP_BUILD_NO_OPTIMIZEGRAPH_PROCESS",
		"ASSIMP_BUILD_NO_SPLITBYBONECOUNT_PROCESS",
		"ASSIMP_BUILD_NO_DEBONE_PROCESS",
		"ASSIMP_BUILD_NO_EMBEDTEXTURES_PROCESS",
		"ASSIMP_BUILD_NO_GLOBALSCALE_PROCESS",
  }

  files 
	{
		-- Dependencies
		'${IncludeDir.assimp}/assimp/contrib/unzip/**',
		'${IncludeDir.assimp}/assimp/contrib/irrXML/**',
		'${IncludeDir.assimp}/assimp/contrib/zlib/*',
		-- Common
		'${IncludeDir.assimp}/assimp/code/Common/**',
		'${IncludeDir.assimp}/assimp/code/PostProcessing/**',
		'${IncludeDir.assimp}/assimp/code/Material/**',
		'${IncludeDir.assimp}/assimp/code/CApi/**',
		-- Importers
		'${IncludeDir.assimp}/assimp/code/Collada/**',
		'${IncludeDir.assimp}/assimp/code/Obj/**',
		-- 'assimp/code/Blender/**', 'assimp/contrib/poly2tri/poly2tri/**',
		'${IncludeDir.assimp}/assimp/code/FBX/**',
		-- 'assimp/code/glTF2/**',
		-- 'assimp/code/glTF/**',
		'${IncludeDir.assimp}/assimp/code/Assbin/**' -- For caching
  }

	-- Same Eigen configuration as PD-Core (see the header).
	forceincludes { "pd_eigen_config.h" }

	links
	{
		"PD-Core",
		"OpenGL-Core",
    "assimp-vc142-mtd"
	}

	filter "system:windows"
		systemversion "latest"

		defines
		{
			"GLCORE_PLATFORM_WINDOWS"
		}
		buildoptions { "/openmp" }

	filter "configurations:Debug"
		defines "GLCORE_DEBUG"
		runtime "Debug"
		symbols "on"

	filter "configurations:Release"
		defines "GLCORE_RELEASE"
		runtime "Release"
        optimize "on"
