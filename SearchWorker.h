#pragma once

#include <thread>
#include <vector>
#include <deque>
#include <mutex>
#include <Windows.h>

template<class T>
class SearchWorker
{
	std::vector<std::thread> workers;
	std::deque<T*> jobs;
	std::mutex jobsQueueMutex;

public:
	SearchWorker(int pollSize);
	~SearchWorker();

	void workerFunction();

	void addJob(T* job);
	bool allJobDone();
};




template<class T>
SearchWorker<T>::SearchWorker(int pollSize)
{
	for (int i = 0; i < pollSize; i++)
	{
		workers.push_back(std::thread(&SearchWorker::workerFunction, this));
	}
}

template<class T>
SearchWorker<T>::~SearchWorker()
{
}

template<class T>
void SearchWorker<T>::workerFunction()
{
	while (1)
	{
		T * job = NULL;

		jobsQueueMutex.lock();
		if (jobs.size() > 0)
		{
			job = jobs.front();
			jobs.pop_front();
		}
		jobsQueueMutex.unlock();

		if (job == NULL)
		{
			Sleep(1);
		}
		else
		{
			job->doJob();
		}
	}
}

template<class T>
void SearchWorker<T>::addJob(T * job)
{
	jobsQueueMutex.lock();
	jobs.push_back(job);
	jobsQueueMutex.unlock();
}

template<class T>
bool SearchWorker<T>::allJobDone()
{
	bool result;
	jobsQueueMutex.lock();
	result = jobs.size() == 0;
	jobsQueueMutex.unlock();
	return result;
}

