#include<vector>
#ifndef FLOWCOMPONENT
#define FLOWCOMPONENT


class volComponent;

class flowLink
{
public:
	double DMDT;

	// 连接件上游和下游只有一个元件
	volComponent* comNext;
	volComponent* comLast;

	flowLink() { comNext = nullptr; comLast = nullptr; };

	//连接两个部件
	virtual bool connect(volComponent* connectFrom, volComponent* connectTo)=0;

	/*virtual double get_DMDT() = 0;*/

};


class volComponent
{
public:
	//可能有多个阀门连接
	std::vector<flowLink*> flowNext;
	std::vector<flowLink*> flowLast;

	
	/*flowLink** flowNext;
	flowLink** flowLast;*/

	//虚函数，更新容积内的状态
	//virtual bool update() = 0;

	//由上游节点或者下游节点获取该节点的焓值或者广义过量空气系数
	virtual double get_h() const = 0;

	virtual double get_AFA() const = 0;

	virtual double get_p() const= 0;

	virtual double get_T() const= 0;

	virtual double get_R() const= 0;

	virtual double get_k() const = 0;

};

#endif // !FLOWCOMPONENT


