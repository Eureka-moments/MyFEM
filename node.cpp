#include "node.h"

Node::Node()
{
    nodenum = 0;
    X = 0;
    Y = 0;
    Fx = 0;
    Fy = 0;

}

void Node::set(ULL num, double x, double y) {
    nodenum = num;
    X = x;
    Y = y;
}