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
	while(this->turns[section].turn != waitFor) ;
}

void order::end(const unsigned int section)
{
	if(!ok) return;
	
	if(this->turns.find(section) == this->turns.end())
	{
		std::cerr << "Must start(" << section << ") before you can end it!" << std::endl;
		throw new std::exception();
	}
	
	#pragma omp critical
	{
		this->turns[section].turn++;
	}
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
