// prodconsA2.sc: Producer-Consumer Example in SpecC
// EECS 222 solution A2

#include <stdio.h>

interface IP
{
  void Send(char);
};

interface IC
{
  char Receive(void);
};

channel C implements IP, IC
{
  event Req;
  char Data;
  event Ack;

  void Send(char X)
  {
    Data = X;
    notify Req;
    wait Ack;
  }

  char Receive(void)
  {
    char Y;
    wait Req;
    Y = Data;
    notify Ack;
    return Y;
  }
};

behavior Prod(IP Port)
{
  const char *Msg = "Apples and Oranges";
  const char *p;
  void main(void)
  {
    printf("Producer starts.\n");
    p = Msg;
    do
    { printf("Producer sends '%c'.\n", *p);
      Port.Send(*p);
      p++;
    } while(*p != 0);
    Port.Send(*p);
    printf("Producer done.\n");
  }
};

behavior Cons(IC Port)
{
  char Y;
  void main(void)
  {
    printf("Consumer starts.\n");
    while(true)
    { Y = Port.Receive();
      if (Y == 0)
	  break;
      printf("Consumer received '%c'.\n", Y);
    }
    printf("Consumer done.\n");
  }
};

behavior Main(void)
{
  C c;
  Prod prod(c);
  Cons cons(c);

  int main(void)
  {
    printf("Main starts.\n");
    par{ prod;
         cons;
        }
    printf("Main done.\n");
    return 0;
  }
};

// EOF

