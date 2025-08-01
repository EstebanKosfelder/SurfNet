
namespace SurfNet
{
    using System;
    using static DebugLog;
    public abstract class BasicFace
    {
        public override string ToString()
        {
            return $"{Id:D3}";
        }

        public void Assert()
        {

           this.Halfedge.Circulation();
        }
        public void FixFaces()
        {
            foreach (var he in Halfedge.Circulation())
            {
                    he.Face = this;
               
               
            }
        }
        

        public int Id { get; set; } 
      
        public BasicHalfedge Halfedge { get;  set; }

        public BasicFace(int id) 
        {
            Id = id;
          
        }


      
    }
}