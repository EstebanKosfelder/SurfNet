namespace SurfNet
{
    public class CollapseEvent : CollapseSpec
    {

        public override string ToString()
        {
            return $"kt:{t.Id} {base.ToString()}";
        }
        public KineticTriangle t;
        public CollapseSpec get_priority() => this;
        public CollapseEvent(KineticTriangle t, double now) : base(t.get_collapse(now))
        {

            this.t = t;

        }

        public void update_collapse(double now)
        {

            if(t.is_collapse_spec_valid())
            {

            }
            this.update_collapse(t.get_collapse(now));


        }
        public override bool Equals(object? obj)
        {
            bool result = base.Equals(obj);

            //if (result && obj is CollapseEvent evnt)
            //{
            //    result = t == evnt.t;
            //}

            return result;

        }
    }
}



