#ifndef TRIPLE_H_INCLUDED
#define TRIPLE_H_INCLUDED


class triple {
public:
    int x;
    int y;
    int z;
    bool operator<(const triple &other) const {
        if (x < other.x) return true;
        if (other.x < x) return false;
        if (y < other.y) return true;
        if (other.y < y) return false;
        return z < other.z;
    }
};




#endif // TRIPLE_H_INCLUDED
