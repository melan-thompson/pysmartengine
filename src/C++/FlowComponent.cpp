#include "FlowComponent.h"

//bool flowLink::connect(volComponent* connectFrom, volComponent* connectTo)
//{
//	//上游的下游节点为本节点
//	connectFrom->flowNext.push_back(*this);
//
//	//下游的上游节点为本节点
//	connectTo->flowLast.push_back(*this);
//
//	return true;
//}

//template<class T>
//bool flowLink::connect(T* connectFrom, T* connectTo)
//{
//	//连接两个对象
//	volComponent* connectFrom1 = (volComponent*)connectFrom;
//
//	volComponent* connectTo1 = (volComponent*)connectTo;
//
//	//上游的下游节点为本节点,连接的上游节点为
//	connectFrom1->flowNext.push_back(*this);
//	comLast = connectFrom1;
//
//	//下游的上游节点为本节点
//	connectTo1->flowLast.push_back(*this);
//	comNext = connectTo1;
//	return true;
//}
