#pragma once

#include <concepts>
#include <iterator>
#include <type_traits>

template<std::integral _T>
struct TriangularIndices
{
	using T = std::make_signed_t<_T>;

	T n;
	TriangularIndices(_T n) : n(n) {}

	struct TriangularIndex
	{
		using Idx = TriangularIndex;
		using value_type = Idx;
		using difference_type = T;

		T i;
		T j;

		T to_n() const
		{
			return (i * (i + 1)) / 2 + j;
		}

		TriangularIndex() = default;
		TriangularIndex(T i, T j) : i{ i }, j{ j } {}

		explicit TriangularIndex(T n)
		{
			i = (sqrt(1 + 8 * n) - 1) / 2;
			while ((i * (i + 1)) / 2 > n) i--;
			while (((i + 1) * (i + 2)) / 2 < n) i++;
			j = n - (i * (i + 1)) / 2;
		}

		Idx operator*() const
		{
			return *this;
		}

		Idx& operator++()
		{
			if (j == i)
			{
				j = 0;
				i++;
			}
			else
			{
				j++;
			}
			return *this;
		}

		Idx operator++(int)
		{
			Idx copy{ *this };
			++(*this);
			return copy;
		}

		Idx& operator--()
		{
			if (j == 0)
			{
				i--;
				j = i;
			}
			else
			{
				j--;
			}
			return *this;
		}

		Idx operator--(int)
		{
			Idx copy{ *this };
			--(*this);
			return copy;
		}

		T operator-(const Idx& rhs) const
		{
			return to_n() - rhs.to_n();
		}

		Idx& operator+=(T n)
		{
			*this = Idx(to_n() + n);
			return *this;
		}

		Idx& operator-=(T n)
		{
			return operator+=(-n);
		}

		Idx operator+(T n) const
		{
			Idx copy{ *this };
			copy += n;
			return copy;
		}

		Idx operator[](T n) const
		{
			return operator+(n);
		}

		Idx operator-(T n) const
		{
			return operator+(-n);
		}

		friend Idx operator+(T n, const Idx& rhs)
		{
			Idx copy{ rhs };
			copy += n;
			return copy;
		}

		auto operator<=>(const Idx&) const = default;
	};

	static_assert(std::random_access_iterator<TriangularIndex>);

	TriangularIndex begin() const
	{
		return { 0, 0 };
	}

	TriangularIndex end() const
	{
		return { n, 0 };
	}
};