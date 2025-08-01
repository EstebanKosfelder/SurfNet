using System;

namespace SurfNetPolyBool
{
	using SurfNet;
	using static SurfNet.Mathex;

    public class LinkedList
    {
        public LinkedList()
        {
            Root = new Node { IsRoot = true };
        }

        private Node Root { get; }

        public bool Exists(Node node)
        {
            if (node == null || Equals(node, Root))
                return false;
            return true;
        }

        public bool IsEmpty()
        {
            return Root.Next == null;
        }

        public Node GetHead()
        {
            return Root.Next;
        }
        
        public void InsertBefore(Node node, Func<Node, bool> check)
        {
            var last = Root;
            var here = Root.Next;
            while (here != null)
            {
                if (check(here))
                {
                    node.LinkPrev(here.Prev);
                    node.LinkNext(here);
                    
                    return;
                }
                last = here;
                here = here.Next;
            }
            last.Next = node;
            node.Prev = last;
            node.Next = null;
        }

        public Transition FindTransition(Func<Node, bool> check)
        {
            var previous = Root;
            var here = Root.Next;
            while (here != null)
            {
                if (check(here))
                    break;
                previous = here;
                here = here.Next;
            }
            return new Transition
            {
                Before = Equals(previous, Root) ? null : previous,
                After = here,
                Insert = node =>
                {
                    node.LinkPrev(previous);
                    node.LinkNext(here);
                    
                    return node;
                }
            };
        }

        public static Node Node(Node data)
        {
            data.Prev = null;
            data.Next = null;
            data.Remove = () =>
            {
                data.Prev.Next = data.Next;
                if (data.Next != null)
                    data.Next.Prev = data.Prev;
                data.Prev = null;
                data.Next = null;
            };
            return data;
        }
    }
}