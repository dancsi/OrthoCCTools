#pragma once

#include <functional>
#include <iterator>

//SUITABLE ONLY FOR STRING KEYS WITH UP TO 3 CHARACTERS
//NOTE THAT THE INTERFACE IS DIFFERENT FROM std::map

template<
	class T,
	size_t N,
	class Hash = std::hash<std::string>
> 

struct fast_map
{
	T storage[N];
	std::string keys[N];
	using key_type = std::string;
	typedef T value_type;

	class iterator
	{
	private:
		std::string key;
		T* ptr;
	public:
		std::pair<std::string, T&> operator*()
		{
			return{ key, *ptr };
		}
		const bool operator!=(const iterator& rhs) const
		{
			return ptr != rhs.ptr;
		}
		iterator(std::string key, T* ptr) : key(key), ptr(ptr) {}
	};

	iterator end()
	{
		return{ "", storage + N };
	}

	iterator find(const std::string& key)
	{
		size_t h = Hash()(key);
		if (h >= N) return end();
		if (keys[h] != key) return end();
		return{key, storage+h};
	}
	void insert(std::pair<std::string, T>& p)
	{
		size_t h = Hash()(p.first);
		if (h >= N) return;
		keys[h] = p.first;
		storage[h] = p.second;
	}
	void insert(std::pair<std::string, T>&& p)
	{
		size_t h = Hash()(p.first);
		if (h >= N) return;
		keys[h] = p.first;
		storage[h] = p.second;
	}
};