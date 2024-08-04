event Metric(i = 0)
{
    if (is_constant(cm))
    {
        scalar *l = list_copy(all);
        cm = new scalar;
        free(all);
        all = list_concat({cm}, l);
        free(l);
    }

if (is_constant(fm.x))
{
    scalar *l = list_copy(all);
    fm = new face vector;
    free(all);
    all = list_concat((scalar *){fm}, l);
    free(l);
}
scalar cmv = cm;
foreach ()
    cmv[] = 0.1;

face vector fmv = fm;
foreach_face()
{
    fmv.x[] = 1.;
    fmv.y[] = .1;
    // fmv.z[] = 1.; For 3D
}
}