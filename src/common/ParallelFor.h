#pragma once

#include <algorithm>
#include <concepts>
#include <cstdlib>
#include <iterator>
#include <thread>
#include <vector>

template<std::random_access_iterator It, std::invocable<std::iter_value_t<It>> F>
void parallel_for(It first, It last, F f)
{
	int nthreads = std::thread::hardware_concurrency();

	auto work_func = [&](It thread_first, It thread_last) 
	{
		std::for_each(thread_first, thread_last, f);
	};

	if (nthreads == 1)
	{
		work_func(first, last);
		return;
	}

	auto n = std::distance(first, last);
	auto work_size = n / (nthreads - 1);

	std::vector<std::pair<It, It>> work_ranges;
	for (decltype(n) i = 0; i < nthreads - 1; i++)
	{
		work_ranges.push_back({ first + i * work_size, first + (i + 1) * work_size });
	}
	work_ranges.push_back({ first + (nthreads - 1) * work_size, last });

	std::vector<std::thread> threads;
	for (const auto& r: work_ranges)
	{
		auto [thread_first, thread_last] = r;
		threads.emplace_back(work_func, thread_first, thread_last);
	}

	for (auto&& t : threads) 
	{
		t.join();
	}
}