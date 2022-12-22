#ifndef BARSIKSPRITE_H
#define BARSIKSPRITE_H

#include <QString>
#include <QColor>
#include <QImage>
#include <QPixmap>
#include <QDebug>

class BarsikSprite
{
public:
    BarsikSprite()
    {

    }

    BarsikSprite(int w, int h)
    {
        Create(w, h);
    }

    BarsikSprite(QString sFile)
    {
        if (!Load(sFile))
            Create(8, 8);
    }

    int nWidth = 0;
    int nHeight = 0;

    wchar_t *m_Glyphs = nullptr;
    QColor *m_Colours = nullptr;

private:
    void Create(int w, int h)
    {
        nWidth = w;
        nHeight = h;
        m_Glyphs = new wchar_t[w*h];
        m_Colours = new QColor[w*h];
        for (int i = 0; i < w*h; i++)
        {
            m_Glyphs[i] = L' ';
            m_Colours[i] = Qt::red;
        }
    }

public:
    void SetGlyph(int x, int y, wchar_t c)
    {
        if (x <0 || x >= nWidth || y < 0 || y >= nHeight)
            return;
        else
            m_Glyphs[y * nWidth + x] = c;
    }

    void SetColour(int x, int y, QColor c)
    {
        if (x <0 || x >= nWidth || y < 0 || y >= nHeight)
            return;
        else
            m_Colours[y * nWidth + x] = c;
    }

    wchar_t GetGlyph(int x, int y)
    {
        if (x <0 || x >= nWidth || y < 0 || y >= nHeight)
            return L' ';
        else
            return m_Glyphs[y * nWidth + x];
    }

    QColor GetColour(int x, int y)
    {
        if (x <0 || x >= nWidth || y < 0 || y >= nHeight)
            return Qt::red;
        else
            return m_Colours[y * nWidth + x];
    }

    wchar_t SampleGlyph(float x, float y)
    {
        int sx = (int)(x * (float)nWidth);
        int sy = (int)(y * (float)nHeight - 1.0f);
        if (sx <0 || sx >= nWidth || sy < 0 || sy >= nHeight)
            return L' ';
        else
            return m_Glyphs[sy * nWidth + sx];
    }

    QColor SampleColour(float x, float y)
    {
        int sx = (int)(x * (float)nWidth);
        int sy = (int)(y * (float)nHeight - 1.0f);
        if (sx <0 || sx >= nWidth || sy < 0 || sy >= nHeight)
            return Qt::red;
        else
            return m_Colours[sy * nWidth + sx];
    }
    bool Load(QString sFile)
    {
        delete[] m_Glyphs;
        delete[] m_Colours;
        nWidth = 0;
        nHeight = 0;

        QImage img;
        img.load(sFile);


qDebug() << " load image" ;
        nWidth = img.width();
        nHeight = img.height();
qDebug() << " creator" ;
        Create(nWidth, nHeight);


qDebug() << " colors array" ;

        for(int x = nWidth; x > 0; x--) {
            for(int y = nHeight; y > 0; y--) {
                int xx = nWidth - x;
                int yy = nHeight - y;
                QColor c(img.pixel(x,y));
                m_Colours[xx + yy*nWidth] = c;
            }
        }

        return true;
    }

};


#endif // BARSIKSPRITE_H
