# Visilib, an open source library for exact visibility computation.
# Copyright(C) 2021 by Denis Haumont
#
# This file is part of Visilib.
#
# Visilib is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Visilib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Visilib. If not, see <http://www.gnu.org/licenses/>

set(GeometrySrc
   geometry_aabbox.h
   geometry_convex_polygon.h
   geometry_ray.h
   geometry_convex_hull.h
   geometry_position_type.h
   geometry_mesh_description.h
   geometry_occluder_set.h
   )

set(HelperSrc
    helper_statistic_collector.h
    helper_geometry_scene_reader.h
    helper_synthetic_mesh_builder.h
    helper_triangle_mesh.h
    helper_triangle_mesh_container.h
	helper_visual_debugger.h
    )

set(MathSrc
    math_arithmetic.h
    math_geometry.h
    math_combinatorial.h
    math_predicates.h
    math_matrix_4.h
    math_plucker_2.h
    math_plucker_6.h
    math_vector_2.h
    math_vector_3.h
    math_plane_3.h
	math_gmp_float.h
	math_gmp_rational.h
	math_real_expr.h
	math_mpfr.h
    )

set(PluckerSrc
    plucker_polyhedron.h
    plucker_polytope.h
    plucker_polytope_builder.h
    plucker_polytope_complex.h
    plucker_polytope_splitter.h
    )

set(VisibilitySrc
    visibility_solver.h
    visibility_aperture_finder.h
    visibility_exact_query.h
	visibility_ray.h
    )

set(SilhouetteSrc
    silhouette_processor.h
    silhouette.h
    silhouette_container.h
    silhouette_container_embree.h
    silhouette_mesh_face.h
    )
set(ConvexHullSrc
    external/convexHull/convexHull.h
    )

set(Headers
    visilib.h
    visilib.hpp
    visilib_core.h
    )

source_group( geometry FILES ${GeometrySrc} )
source_group( helper FILES ${HelperSrc} )
source_group( math FILES ${MathSrc} )
source_group( plucker FILES ${PluckerSrc} )
source_group( silhouette FILES ${SilhouetteSrc} )
source_group( visibility FILES ${VisibilitySrc} )
source_group( external FILES ${ConvexHullSrc} )
source_group( include FILES ${Headers} )
source_group( src FILES ${Sources} )

add_library( ${VISILIB_LIB_NAME} STATIC
    ${GeometrySrc}
    ${GeometryHeader}
    ${HelperSrc}
    ${MathSrc}
    ${PluckerSrc}
    ${SilhouetteSrc}
    ${VisibilitySrc}
    ${ConvexHullSrc}
    ${RootSrc}
    ${Headers}
    )

set_target_properties(${VISILIB_LIB_NAME} PROPERTIES LINKER_LANGUAGE CXX)
