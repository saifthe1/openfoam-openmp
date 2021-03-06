#include "timestamp.hpp"

ts* ts::instance = NULL;

void ts::start(const std::string& stamp_name)
{
	if(!ok) return;
	
	if(this->stamps[stamp_name].flip == false)
	{
		std::cerr << "Can't start(" << stamp_name << "), since the last one hasn't ended yet!" << std::endl;
		throw new std::exception();
	}
	this->stamps[stamp_name].flip = false;
	
	if(this->stamps.find(stamp_name) == this->stamps.end()) // Didn't find stamp_name
	{
		this->stamps[stamp_name] = ts_obj(::omp_get_wtime());
	}
	else
	{
		this->stamps[stamp_name].lastStamp = ::omp_get_wtime();
	}
}

void ts::end(const std::string& stamp_name)
{
	if(!ok) return;

	if(this->stamps.find(stamp_name) == this->stamps.end())
	{
		std::cerr << "Must start(" << stamp_name << ") before you can end it!" << std::endl;
		throw new std::exception();
	}
	
	if(this->stamps[stamp_name].flip == true)
	{
		std::cerr << "Can't end(" << stamp_name << "), since it's not started!" << std::endl;
		throw new std::exception();
	}
	this->stamps[stamp_name].flip = true;
	
	this->stamps[stamp_name].time += omp_get_wtime() - this->stamps[stamp_name].lastStamp;
	this->stamps[stamp_name].stamps += 1;
}

void ts::toggle()
{
	this->ok = !this->ok;
}
void ts::toggle(bool set)
{
	this->ok = set;
}

const std::map<std::string, ts_obj> ts::getStamps()
{
	return stamps;
}

void ts::printStamps()
{
	if(!ok) return;

	std::cout 	<< std::setiosflags(std::ios::left)
				<< std::setw(50) << "What" << " "
				<< std::setw(18) << "Time" << " "
				<< std::setw(8)  << "Count" << std::endl;
	for(std::map<std::string, ts_obj>::iterator iter = stamps.begin(); iter != stamps.end(); ++iter)
	{
		std::cout 	<< std::setiosflags(std::ios::left)
					<< std::setw(50) << std::setfill('.') << iter->first << " " << std::setfill(' ')
					<< std::setw(18) << std::setprecision(5) << iter->second.time << " "
					<< std::setw(8)  << iter->second.stamps << std::endl;
	}
}

ts* ts::get_instance()
{
	return instance ? instance : (instance = new ts);
}
