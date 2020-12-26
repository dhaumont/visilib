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

#include <string>
#include <iostream>

#include "math_combinatorial.h"

using namespace visilib;

bool MathCombinatorialTest(std::string& )
{
	{//MathCombinatorial::hasFacet
		if (MathCombinatorial::hasFacet(std::vector<size_t>(), 1))
			return false;

		if (MathCombinatorial::hasFacet(std::vector<size_t>{ 1, 2, 3, 4 }, 0))
			return false;

		for (size_t i = 1; i <= 6; i++)
		{
			if (!MathCombinatorial::hasFacet(std::vector<size_t>{ 1, 2, 3, 4, 5, 6 }, i))
				return false;
		}
	}

	{//MathCombinatorial::haveAtLeastNCommonFacets
		for (size_t i = 0; i < 5; i++)
		{
			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 2, 3}, std::vector<size_t>{1, 2, 3}, i) != (i <= 3))
				return false;

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 2, 3, 4, 5, 6}, std::vector<size_t>{1, 2, 3}, i) != (i <= 3))
				return false;

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 2, 3}, std::vector<size_t>{4, 5, 6, 7, 8, 9 }, i) != (i <= 0))
				return false;

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 6, 9, 12, 15, 24, 31}, std::vector<size_t>{6, 8, 10, 22, 24, 25, 31, 32, 33 }, i) != (i <= 3))
				return false;

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 12, 13, 14, 15, 16}, std::vector<size_t>{1, 2, 3, 4, 5, 6, 7, 16  }, i) != (i <= 2))
				return false;

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{5, 10, 15, 20, 25, 30}, std::vector<size_t>{11, 12, 13, 14, 15, 16, 17, 18,19  }, i) != (i <= 1))
				return false;

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 1, 1, 1}, std::vector<size_t>{1, 1, 1, 1, 1, 1 }, i) != (i <= 4))
				return false;

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{}, std::vector<size_t>{1, 2, 3, 4 }, i) != (i <= 0))
				return false;

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 2, 3, 4}, std::vector<size_t>{}, i) != (i <= 0))
				return false;

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{}, std::vector<size_t>{}, i) != (i <= 0))
				return false;
		}
	}

	{//MathCombinatorial::initFacets
		std::vector<size_t> result;
		MathCombinatorial::initFacets(std::vector<size_t>{1, 2, 3}, std::vector<size_t>{1, 2, 3}, 4, result);
		if (result != std::vector<size_t>{ 1, 2, 3, 4 })
			return false;

		MathCombinatorial::initFacets(std::vector<size_t>{1, 2, 3}, std::vector<size_t>{1, 2, 3, 4, 5, 6, 7}, 8, result);
		if (result != std::vector<size_t>{ 1, 2, 3, 8 })
			return false;

		MathCombinatorial::initFacets(std::vector<size_t>{1, 2, 3, 4, 5, 6, 7}, std::vector<size_t>{1, 2, 3}, 8, result);
		if (result != std::vector<size_t>{ 1, 2, 3, 8 })
			return false;

		MathCombinatorial::initFacets(std::vector<size_t>{1, 2, 3, 4, 5, 6}, std::vector<size_t>{2, 4, 6, 8, 10}, 11, result);
		if (result != std::vector<size_t>{ 2, 4, 6, 11 })
			return false;

		MathCombinatorial::initFacets(std::vector<size_t>{1, 2, 3, 4, 5, 6}, std::vector<size_t>{1, 2, 3, 4, 5, 6}, 7, result);
		if (result != std::vector<size_t>{ 1, 2, 3, 4, 5, 6, 7 })
			return false;
	}

	return true;
}
bool MathPluckerTest(std::string& )
{
	std::cout << "Hello" << std::endl;
	//    areVisible(nullptr, nullptr, 0, nullptr, 0);

	return true;
}

int main()
{
	std::string errorMessage;

	if (!MathCombinatorialTest(errorMessage))
	{
		return 1;
	}

	return 0;
}