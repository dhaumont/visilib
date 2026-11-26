/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2020 by Denis Haumont

This file is part of Visilib.

Visilib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Visilib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Visilib. If not, see <http://www.gnu.org/licenses/>
*/


bool MathCombinatorialTest(std::string&);
bool MathPluckerParallelTest(std::string&);
bool MathBeamCapsuleTest(std::string&);
bool MathShadowGeometryTest(std::string&);
bool MathIncidenceTest(std::string&);
bool MathIntegerPredicateTest(std::string&);
bool MathNormalizationTest(std::string&);
bool testConfiguration(bool& retflag);
bool VisibilityTest(std::string&);
bool DirectionalVisibilityTest(std::string&);
bool SolverConsistencyTest(std::string&);
bool SolverFatRayTest(std::string&);
bool SolverEdgeCases(std::string&);
bool SolverDegenerateTests(std::string&);
bool SolverModelTests(std::string&);
bool SolverPrecisionTests(std::string&);
bool ShadowDiscontinuityTest(std::string&);

void selectSolverTestVariants(const std::string& idList);

const double cDefaultQueryBudgetSeconds = 2.0;

/* The models of demo/data are read whole, and the big ones are answered by a search that has to
walk a polytope complex of tens of thousands of triangles before it can say anything. The two
seconds the small synthetic scenes are graded with stops them before the scene is even built, so
their suite asks for ten. */
const double cModelQueryBudgetSeconds = 10.0;

double queryBudgetSeconds();
void setQueryBudgetSeconds(double seconds);
void disableQueryBudget();
