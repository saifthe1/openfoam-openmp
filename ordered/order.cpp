#include "order.hpp"

order* order::instance = NULL;

void order::start(const unsigned int section, const unsigned int waitFor)
{
	if(!ok) return;
	
	if(this->turns.find(section) == this->turns.end()) // Didn't find section
	{
		// Add an order object
		this->turns[section] = order_obj();
	}

	#pragma omp critical
	std::cout << "Thread #" << omp_get_thread_num() << " Queues Sect#" << section << " [" << this->turns[section].turn << "/" << waitFor << "]" << std::endl;

	while(this->turns[section].turn != waitFor) ;

	#pragma omp critical
	std::cout << "Thread #" << omp_get_thread_num() << " Enters Sect#" << section << std::endl;
}

void order::end(const unsigned int section)
{
	if(!ok) return;
	
	#pragma omp critical
	std::cout << "Thread #" << omp_get_thread_num() << " Leaves Sect#" << section << std::endl;
	
	if(this->turns.find(section) == this->turns.end())
	{
		std::cerr << "Must start(" << section << ") before you can end it!" << std::endl;
		throw new std::exception();
	}
	
	#pragma omp critical
	{
		this->turns[section].turn++;
	}
	std::cout << "Thread #" << omp_get_thread_num() << " Leaves Sect#" << section << " with turn=" << this->turns[section].turn << std::endl;
}

void order::toggle()
{
	this->ok = !this->ok;
}
void order::toggle(bool set)
{
	this->ok = set;
}

order* order::get_instance()
{
	return instance ? instance : (instance = new order);
}
