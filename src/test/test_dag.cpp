#include <DagAligner.hpp>

int main(){
    vector<xy_match> matches = {
            {{0,0}, 1},
            {{0,1}, 1},
            {{2,2}, 1},
            {{2,3}, 1},
            {{3,4}, 1},
    };

    Dag dag(matches,1);

    return 0;
}
