
// Includes
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <map>

#define ORDER_START order::get_instance()->start
#define ORDER_END order::get_instance()->end
#define ORDER_TOGGLE order::get_instance()->toggle
#define ORDER_PRINT order::get_instance()->printStamps

class order_obj
{
	public:
		unsigned int turn;
		bool flip;
		order_obj() : turn(0), flip(true) {}
};

class order
{
	private:
		
		std::map<unsigned int, order_obj> turns;
		
		// Singelton
		order() : ok(true) {};
		order(const order&) : ok(true) {};
		static order* instance;

		bool ok;
		
	public:
		void start(const unsigned int section, const unsigned int waitFor);
		void end(const unsigned int section);

		void toggle();
		void toggle(bool);

		static order* get_instance();
};


