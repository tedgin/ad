import std.stdio;
import ad.math;
import std.traits;

alias GDN = GenDualNum;

struct S(ulong Degree)
{
    enum ulong DEGREE = Degree;

    int x;

    this(int v) { x = v; }

    T opCast(T)() const if (fullyQualifiedName!(TemplateOf!T) == "main.S")
    {
        return S!(T.DEGREE)(1);
    }
}

void main()
{
    auto x = GenDualNum!1(2);
    auto y = cast(GenDualNum!2) x;
    assert(y == x && typeof(y).DEGREE == 2);
}