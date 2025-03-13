/*******************************************************************************
 * Copyright (c) 2016, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.common.networksimplex;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

/**
 * A node used by the {@link NetworkSimplex} algorithm. It has a set of incoming and outgoing edges.
 */
public final class NNode {

    // SUPPRESS CHECKSTYLE NEXT 15 VisibilityModifier
    /** A public id, unused internally, use it for whatever you want. */
    public int id;
    /** Internally set and used id to index arrays. */
    protected int internalId;
    /** An object from which this edge is derived. */
    public Object origin;
    /** The type is attached as label to the debug graph. Apart from this it has no semantic meaning. */
    public String type = "";

    /**
     * The layer each node is currently assigned to. Note that during layerer execution, the lowest
     * layer is not necessary the zeroth layer. To fulfill this characteristic, a final call of
     * {@link NetworkSimplex#normalize()} has to be performed.
     */
    public int layer;
    
    // In addition to a node's incoming and outgoing edges we cache a list of 
    // the union of these two. To be able to tell when we have to update this list,
    // we listen (sort of) to changes to the two underlying lists. 
    private ChangeAwareArrayList<NEdge> outgoingEdges = new ChangeAwareArrayList<>();
    private int outgoingEdgesModCnt = -1;
    private ChangeAwareArrayList<NEdge> incomingEdges = new ChangeAwareArrayList<>();
    private int incomingEdgesModCnt = -1;
    
    /** Internally cached list of all edges. */
    private ArrayList<NEdge> allEdges = Lists.newArrayList();
    
    
    // SUPPRESS CHECKSTYLE NEXT 20 VisibilityModifier 
    /**
     * A flag indicating whether a specified node is part of the spanning tree determined by
     * {@code tightTree()}.
     * 
     * @see NetworkSimplex#tightTreeDFS() 
     */
    protected boolean treeNode = false;

    /**
     * A collection of edges for which cutvalues are unknown. Cutvalues are updated during every
     * iteration of the network simplex algorithm.
     * 
     * In principle, no order is required. Still, we use {@link ArrayList} because experiments
     * showed, that an {@link ArrayList} performs faster here than a {@link java.util.HashSet
     * HashSet}.
     */
    protected Collection<NEdge> unknownCutvalues = Lists.newArrayList();

    private NNode() { }
    
    /**
     * @return a {@link NNodeBuilder} to construct a new node.
     */
    public static NNodeBuilder of() {
        return new NNodeBuilder();
    }
    
    /**
     * @return the outgoingEdges
     */
    public List<NEdge> getOutgoingEdges() {
        return outgoingEdges;
    }
    
    /**
     * @return the incomingEdges
     */
    public List<NEdge> getIncomingEdges() {
        return incomingEdges;
    }
    
    /**
     * @return a list with the union of {@link #getOutgoingEdges()} and {@link #getIncomingEdges()}.
     *         The list is cached internally, subsequent calls return the same list as long as
     *         neither of the incoming and outgoing edges changes.
     */
    public List<NEdge> getConnectedEdges() {
        if (incomingEdgesModCnt != incomingEdges.getModCnt()
                || outgoingEdgesModCnt != outgoingEdges.getModCnt()) {
            
            allEdges.clear();
            allEdges.ensureCapacity(incomingEdges.size() + outgoingEdges.size());
            allEdges.addAll(incomingEdges);
            allEdges.addAll(outgoingEdges);

            incomingEdgesModCnt = incomingEdges.getModCnt();
            outgoingEdgesModCnt = outgoingEdges.getModCnt();
        }

        return allEdges;
    }
    
    @Override
    public String toString() {
        return this.origin.toString();
    }
    
    /**
     * Builder class for {@link NNode}s. 
     */
    public static final class NNodeBuilder {
        
        /** The node currently being constructed. */  
        private NNode node;
        
        private NNodeBuilder() {
            node = new NNode();
        }

        /**
         * Sets the id field.
         * 
         * @param id
         *            a non-negative integer.
         * @return this builder.
         */
        public NNodeBuilder id(final int id) {
            node.id = id;
            return this;
        }

        /**
         * Sets the origin field.
         * 
         * @param origin
         *            any object.
         * @return this builder.
         */
        public NNodeBuilder origin(final Object origin) {
            node.origin = origin;
            return this;
        }

        /**
         * Sets the type field.
         * 
         * @param type
         *            any string, has no semantic meaning, can be used for debugging.
         * @return this builder.
         */
        public NNodeBuilder type(final String type) {
            node.type = type;
            return this;
        }

        /**
         * Finally creates this node. That is, the node is added to the passed {@link NGraph}.
         * 
         * @param graph
         *            the {@link NGraph} this node belongs to.
         * @return the created {@link NNode}.
         */
        public NNode create(final NGraph graph) {
            graph.nodes.add(node);
            return node;
        }

    }
    
    /**
     * Delegates all its calls to an internal {@link ArrayList} An easy implementation would just expose the
     * {@link java.util.AbstractList AbstractList}'s {@code modCount} field. However this is not compatible with GWT
     * since GWT's list implementations do not increment the 'modCount' variable.
     * 
     * @param <E>
     *            type of the elements in this collection.
     */
    private static class ChangeAwareArrayList<E> implements List<E> {

        private final ArrayList<E> list = Lists.newArrayList(); 
        private int modCount;
        
        // CHECKSTYLEOFF Parameters
        private int getModCnt() {
            return modCount;
        }

        @Override
        public int size() {
            return list.size();
        }

        @Override
        public boolean isEmpty() {
            return list.isEmpty();
        }

        @Override
        public boolean contains(Object o) {
            return list.contains(o);
        }

        @Override
        public Iterator<E> iterator() {
            return Iterators.unmodifiableIterator(list.iterator());
        }

        @Override
        public Object[] toArray() {
            return list.toArray();
        }

        @Override
        public <T> T[] toArray(T[] a) {
            return list.toArray(a);
        }

        @Override
        public boolean add(E e) {
            modCount++;
            return list.add(e);
        }

        @Override
        public boolean remove(Object o) {
            modCount++;
            return list.remove(o);
        }

        @Override
        public boolean containsAll(Collection<?> c) {
            return list.containsAll(c);
        }

        @Override
        public boolean addAll(Collection<? extends E> c) {
            modCount++;
            return list.addAll(c);
        }

        @Override
        public boolean addAll(int index, Collection<? extends E> c) {
            modCount++;
            return list.addAll(index, c);
        }

        @Override
        public boolean removeAll(Collection<?> c) {
            boolean changed = list.removeAll(c);
            if (changed) {
                modCount++;
            }
            return changed;
        }

        @Override
        public boolean retainAll(Collection<?> c) {
            boolean changed = list.retainAll(c);
            if (changed) {
                modCount++;
            }
            return changed;
        }

        @Override
        public void clear() {
            modCount++;
            list.clear();
        }

        @Override
        public E get(int index) {
            return list.get(index);
        }

        @Override
        public E set(int index, E element) {
            modCount++;
            return list.set(index, element);
        }

        @Override
        public void add(int index, E element) {
            modCount++;
            list.add(index, element);
        }

        @Override
        public E remove(int index) {
            modCount++;
            return list.remove(index);
        }

        @Override
        public int indexOf(Object o) {
            return list.indexOf(o);
        }

        @Override
        public int lastIndexOf(Object o) {
            return list.lastIndexOf(o);
        }

        @Override
        public ListIterator<E> listIterator() {
            // if implemented, must guarantee not to alter the list or to increment modCount
            throw new UnsupportedOperationException();
        }

        @Override
        public ListIterator<E> listIterator(int index) {
            // if implemented, must guarantee not to alter the list or to increment modCount
            throw new UnsupportedOperationException();
        }

        @Override
        public List<E> subList(int fromIndex, int toIndex) {
            return list.subList(fromIndex, toIndex);
        }
        // CHECKSTYLEON Parameters
    }
    
}
