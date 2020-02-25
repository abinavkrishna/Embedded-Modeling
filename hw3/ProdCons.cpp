#include <systemc.h>

class Prod_interface : public sc_interface {
  public:
    virtual void Send(char) = 0;
};

class Cons_interface : public sc_interface {
  public:
    virtual char Receive() = 0;
};

class C : public sc_channel,public Prod_interface, public Cons_interface {
  public:
  C(sc_module_name name) :sc_channel(name){}

  sc_event write_event, read_event;
  char Data;

  void Send(char X)
  {
    Data = X;
    write_event.notify(SC_ZERO_TIME);
    wait(read_event);
  }

  char Receive(void)
  {
    char Y;
    wait(write_event);
    Y = Data;
    read_event.notify(SC_ZERO_TIME);
    return Y;
  }
};

class Prod : public sc_module
{
  public:
  sc_port<Prod_interface> out;
  SC_CTOR(Prod) // the module constructor
  {
    SC_THREAD(main); // start the producer process
  }
  void main(void)
  {
    char Msg[] = "Apples and Oranges";
    char *p;
    std :: cout<<"Producer starts."<<std :: endl;
    p = Msg;
    do
    { std :: cout<<"Producer sends '"<<*p<<"'."<<std :: endl;
      out->Send(*p);
      p++;
    } while(*p != 0);
    out->Send(*p);
    std :: cout<<"Producer done."<<std :: endl;
  }
};

class Cons : public sc_module {
  public:
  char Y;
  sc_port<Cons_interface> in;
  SC_CTOR(Cons) // the module constructor
  {
    SC_THREAD(main); // start the consumer process
  }
  void main(void)
  {
    std :: cout<<"Consumer starts."<<std :: endl;
    while(true)
    { Y = in->Receive();
      if (Y == 0)
	  break;
      std :: cout<<"Consumer received '"<<Y<<"'."<<std :: endl;
    }
    std :: cout<<"Consumer done."<<std :: endl;
  }
};

class top : sc_module {
  public:
  C* c_inst;
  Prod* prod_inst;
  Cons* cons_inst;
  top (sc_module_name name) : sc_module(name){
    c_inst = new C("C1");
    prod_inst = new Prod("Producer1");
    prod_inst->out(*c_inst);
    cons_inst = new Cons("Consumer1");
    cons_inst->in(*c_inst);
  }
};


int sc_main(int argc, char* argv[]) {
  top top1("Top1");
  std :: cout<<"sc_main starts"<<std :: endl;
  sc_start();
  std :: cout<<"sc_main ends"<<std :: endl;
  return 0;
};

// EOF

