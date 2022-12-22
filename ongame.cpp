#include "widget.h"

void Widget::OnGame()   {

    OnTouch();
        // main app logic in timer period 15ms (edit in defines.h)

    for(int i = 0; i < BUTTON_COUNT; i++)    {
        Game.Key[i] = myKey[i];
    }

    Game.update();

    repaint();
}
