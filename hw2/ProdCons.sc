#include<stdio.h>
#include<stdlib.h>

interface Sender {
	void Send(char);
};

interface Receiver {
	char Receive(void);
};

channel C1 implements Sender, Receiver {
	event Req;
	char Data;
	event Ack;

	void Send(char ch) {
	    Data = ch;
	    notify Req;
	    wait Ack;
	}

	char Receive(void) {
		char ch;
		wait Req;
		ch = Data;
		notify Ack;
		return ch;
	} 
};

behavior Producer(Sender Port) {
	char *ch = "Apples and Oranges";
	void main(void) {
		int i = 0;
		printf("Producer starts.\n");
		do {
			printf("Producer sends '%c'. \n",ch[i]);
			Port.Send(ch[i]);
			i++;
		}while(ch[i] != '\0');
		printf("Producer done. \n");
	}
};

behavior Consumer(Receiver Port) {
        char ch;
	void main(void) {
		printf("Consumer starts.\n");
		while(true){
			ch = Port.Receive();
			if(ch == 0) break;
			printf("Consumer received '%c'. \n",ch);
		}
		printf("Consumer done. \n");
	}	
};

behavior Main() {
	C1 c;
	Producer prod(c);
	Consumer cons(c);
	int main() {
		printf("Main starts.\n");
		par{
			prod;
			cons;
		}
		printf("Main done.\n");
		return 0;
	}
};
