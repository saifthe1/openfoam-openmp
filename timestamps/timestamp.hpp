
// Includes
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <map>

#define TS_START ts::get_instance()->start
#define TS_END ts::get_instance()->end
#define TS_TOGGLE ts::get_instance()->toggle
#define TS_PRINT ts::get_instance()->printStamps

class ts_obj
{
	public:
		unsigned int stamps;
		double time, lastStamp;
		bool flip;
		ts_obj() : stamps(0), time(0), lastStamp(0), flip(true) {}
		ts_obj(const double curTime) : stamps(0), time(0), lastStamp(curTime), flip(true) {}
};

class ts
{
	private:
		
		std::map<std::string, ts_obj> stamps;
		
		// Singelton
		ts() : ok(true) {};
		ts(const ts&) : ok(true) {};
		static ts* instance;

		bool ok;
		
	public:
		void start(const std::string& stamp_name);
		void end(const std::string& stamp_name);

		void toggle();
		void toggle(bool);

		const std::map<std::string, ts_obj> getStamps();
		void printStamps();

		static ts* get_instance();
};


