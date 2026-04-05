// ============================================================================
//  main.cpp  —  SolidStateCore entry point
// ============================================================================

#include "GuiApp.h"
#include <iostream>

int main() {
    GuiApp app;
    if (!app.init(1600, 900, "SolidStateCore  —  Condensed Matter Physics Simulator")) {
        std::cerr << "Failed to initialize GUI.\n";
        return 1;
    }
    app.run();
    return 0;
}
