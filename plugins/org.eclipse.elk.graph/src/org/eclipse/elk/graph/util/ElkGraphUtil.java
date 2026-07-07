/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.util;

import java.lang.annotation.Annotation;
import java.util.Arrays;
import java.util.ConcurrentModificationException;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.NoSuchElementException;
import java.util.Objects;

import org.eclipse.elk.graph.EMapPropertyHolder;
import org.eclipse.elk.graph.ElkBendPoint;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkGraphFactory;
import org.eclipse.elk.graph.ElkGraphPackage;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.properties.AdvancedPropertyValue;
import org.eclipse.elk.graph.properties.ExperimentalPropertyValue;
import org.eclipse.emf.common.util.AbstractTreeIterator;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.EStructuralFeature;
import org.eclipse.emf.ecore.util.EContentsEList.FeatureFilter;
import org.eclipse.emf.ecore.util.EContentsEList.FeatureIteratorImpl;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Utility methods that make using the ElkGraph data structure a bit easier and thereby make ELK great again! The
 * class offers different types of methods described below.
 *
 *
 * <h2>Factory Methods</h2>
 *
 * <p>While {@link ElkGraphFactory} offers methods that simply create new model objects, the factory methods
 * offered by this class are designed to make creating a graph programmatically as easy as possible. This
 * includes automatically setting containments and may at some point also include applying sensible defaults,
 * where applicable. The latter is why using these methods instead of factory methods is usually a good
 * idea.</p>
 *
 *
 * <h2>Convenience Methods</h2>
 *
 * <p>This class offers quite a few convenience methods, including easy iteration over end points of an edge
 * or edges incident to a node, a way to find the node a connectable shape represents, as well as ways to
 * interact with a graph's structure.</p>
 *
 * TODO: More documentation about what's in here.
 */
public final class ElkGraphUtil {

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    /**
     * Creates a new root node that will represent a graph.
     *
     * @return the root node.
     */
    public static ElkNode createGraph() {
        return createNode(null);
    }

    /**
     * Creates a new node in the graph represented by the given parent node.
     *
     * @param parent the parent node. May be {@code null}, in which case the new node is not added to anything.
     * @return the new node.
     */
    public static ElkNode createNode(final ElkNode parent) {
        ElkNode node = ElkGraphFactory.eINSTANCE.createElkNode();

        if (parent != null) {
            node.setParent(parent);
        }

        return node;
    }

    /**
     * Creates a new port for the given parent node.
     *
     * @param parent the parent node. May be {@code null}, in which case the new port is not added to anything.
     * @return the new port.
     */
    public static ElkPort createPort(final ElkNode parent) {
        ElkPort port = ElkGraphFactory.eINSTANCE.createElkPort();

        if (parent != null) {
            port.setParent(parent);
        }

        return port;
    }

    /**
     * Creates a new label for the given graph element.
     *
     * @param parent the parent element. May be {@code null}, in which case the new label is not added to anything.
     * @return the new label.
     */
    public static ElkLabel createLabel(final ElkGraphElement parent) {
        ElkLabel label = ElkGraphFactory.eINSTANCE.createElkLabel();

        if (parent != null) {
            label.setParent(parent);
        }

        return label;
    }

    /**
     * Creates a new label with the given text for the given graph element.
     *
     * @param text the label's text.
     * @param parent the parent element. May be {@code null}, in which case the new label is not added to anything.
     * @return the new label.
     */
    public static ElkLabel createLabel(final String text, final ElkGraphElement parent) {
        ElkLabel label = createLabel(parent);
        label.setText(text);
        return label;
    }

    /**
     * Creates a new edge contained in the given node, but not connecting anything yet. Note that the containing
     * node defines the coordinate system for the edge's routes and is not straightforward to select. One way to get
     * around this issue is to create an edge without a containing node, setup its sources and targets, and call
     * {@link #updateContainment(ElkEdge)} afterwards.
     *
     * @param containingNode the edge's containing node. May be {@code null}, in which case the new edge is not added
     *                       to anything.
     * @return the new edge.
     */
    public static ElkEdge createEdge(final ElkNode containingNode) {
        ElkEdge edge = ElkGraphFactory.eINSTANCE.createElkEdge();

        if (containingNode != null) {
            edge.setContainingNode(containingNode);
        }

        return edge;
    }

    /**
     * Creates an edge that connects the given source to the given target and sets the containing node accordingly.
     * This requires the source and target to be in the same graph model.
     *
     * @param source the edge's source.
     * @param target the edge's target.
     * @return the new edge.
     * @throws NullPointerException if {@code source} or {@code target} is {@code null}.
     */
    public static ElkEdge createSimpleEdge(final ElkConnectableShape source, final ElkConnectableShape target) {
        Objects.requireNonNull(source, "source cannot be null");
        Objects.requireNonNull(target, "target cannot be null");

        ElkEdge edge = createEdge(null);

        edge.getSources().add(source);
        edge.getTargets().add(target);
        updateContainment(edge);

        return edge;
    }

    /**
     * Creates a hyperedge that connects the given sources to the given targets and sets the containing node
     * accordingly. This requires the sources and targets to be in the same graph model.
     *
     * @param sources the edge's sources.
     * @param targets the edge's targets.
     * @return the new edge.
     * @throws NullPointerException if {@code sources} or {@code targets} is {@code null}.
     */
    public static ElkEdge createHyperedge(final Iterable<ElkConnectableShape> sources,
            final Iterable<ElkConnectableShape> targets) {

        Objects.requireNonNull(sources, "sources cannot be null");
        Objects.requireNonNull(targets, "targets cannot be null");

        ElkEdge edge = createEdge(null);

        List<ElkConnectableShape> edgeSources = edge.getSources();
        for (ElkConnectableShape source : sources) {
            edgeSources.add(source);
        }

        List<ElkConnectableShape> edgeTargets = edge.getTargets();
        for (ElkConnectableShape target : targets) {
            edgeTargets.add(target);
        }

        updateContainment(edge);

        return edge;
    }

    /**
     * Creates a new edge section and adds it to the given edge.
     *
     * @param edge the edge to add the new section to. May be {@code null}, in which case the new edge section is
     *             not added to anything.
     * @return the new edge section.
     */
    public static ElkEdgeSection createEdgeSection(final ElkEdge edge) {
        ElkEdgeSection section = ElkGraphFactory.eINSTANCE.createElkEdgeSection();

        if (edge != null) {
            edge.getSections().add(section);
        }

        return section;
    }

    /**
     * Returns the edge's first edge section or creates one if there is none.
     *
     * @param edge the edge whose first edge section to retrieve.
     * @param resetSection {@code true} if all bend points should be removed and the location reset on the section
     *                     before it is returned.
     * @param removeOtherSections {@code true} if all the other edge sections, if any, should be removed.
     * @return the first edge section.
     */
    public static ElkEdgeSection firstEdgeSection(final ElkEdge edge, final boolean resetSection,
            final boolean removeOtherSections) {

        if (edge.getSections().isEmpty()) {
            // Create and return a new section
            return createEdgeSection(edge);
        } else {
            // Retrieve the first section
            ElkEdgeSection section = edge.getSections().get(0);

            if (resetSection) {
                section.getBendPoints().clear();
                section.setStartLocation(0, 0);
                section.setEndLocation(0, 0);
            }

            if (removeOtherSections) {
                List<ElkEdgeSection> sections = edge.getSections();
                while (sections.size() > 1) {
                    sections.remove(sections.size() - 1);
                }
            }

            return section;
        }
    }

    /**
     * Creates a new bend point and adds it to the given edge section.
     *
     * @param edgeSection the edge section to add the bend point to. May be {@code null}, in which case the new
     *                    bend point is not added to anything.
     * @return the new bend point.
     */
    public static ElkBendPoint createBendPoint(final ElkEdgeSection edgeSection) {
        return createBendPoint(edgeSection, 0, 0);
    }

    /**
     * Creates a new bend point with the given coordinates and adds it to the given edge section.
     *
     * @param edgeSection the edge section to add the bend point to. May be {@code null}, in which case the new
     *                    bend point is not added to anything.
     * @param x the bend point's x coordinate.
     * @param y the bend point's y coordinate.
     * @return the new bend point.
     */
    public static ElkBendPoint createBendPoint(final ElkEdgeSection edgeSection, final double x, final double y) {
        ElkBendPoint bendPoint = ElkGraphFactory.eINSTANCE.createElkBendPoint();

        bendPoint.set(x, y);

        if (edgeSection != null) {
            edgeSection.getBendPoints().add(bendPoint);
        }

        return bendPoint;
    }



    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Edge Containment

    /**
     * Changes an edge's containment to the one returned by {@link #findBestEdgeContainment(ElkEdge)}.
     *
     * @param edge the edge to update the containment of.
     * @throws NullPointerException if {@code edge} is {@code null}.
     * @throws IllegalArgumentException if {@code edge} does not have at least one source or target.
     */
    public static void updateContainment(final ElkEdge edge) {
        Objects.requireNonNull(edge, "edge cannot be null");

        edge.setContainingNode(findBestEdgeContainment(edge));
    }

    /**
     * Finds the node the given edge should best be contained in given the edge's sources and targets. This is usually
     * the first common ancestor of all sources and targets. Finding this containment requires all sources and targets
     * to be contained in the same graph model. If that is not the case, this method will return {@code null}. If the
     * edge is not connected to anything, an exception is thrown. If the edge is connected to only one shape, that
     * shape's parent is returned.
     *
     * @param edge the edge to find the best containment for.
     * @return the best containing node, or {@code null} if none could be found.
     * @throws NullPointerException if {@code edge} is {@code null}.
     * @throws IllegalArgumentException if {@code edge} does not have at least one source or target.
     */
    public static ElkNode findBestEdgeContainment(final ElkEdge edge) {
        Objects.requireNonNull(edge, "edge cannot be null");

        /* We start with corner cases: an edge which is not connected to anything and an edge which is connected to
         * just one shape.
         */

        switch (edge.getSources().size() + edge.getTargets().size()) {
        case 0:
            // We need at least one to work with
            throw new IllegalArgumentException("The edge must have at least one source or target.");

        case 1:
            // Return the parent of the only incident thing
            if (edge.getSources().isEmpty()) {
                return connectableShapeToNode(edge.getTargets().get(0)).getParent();
            } else {
                return connectableShapeToNode(edge.getSources().get(0)).getParent();
            }
        }

        /* From now on we know that the edge has at least two incident shapes. We now check the simple common cases
         * first before moving on to the more complex general case. The most common cases are that the edge is not a
         * hyperedge and one of the following is true:
         *   1. The edge's source and target have the same parent node. In this case, the containment is that parent
         *      node.
         *   2. The source is the target's parent (or the other way around). In that case, the source or the target
         *      is the containment, respectively.
         */

        if (edge.getSources().size() == 1 && edge.getTargets().size() == 1) {
            ElkNode sourceNode = connectableShapeToNode(edge.getSources().get(0));
            ElkNode targetNode = connectableShapeToNode(edge.getTargets().get(0));

            if (sourceNode.getParent() == targetNode.getParent()) {
                return sourceNode.getParent();
            } else if (sourceNode == targetNode.getParent()) {
                return sourceNode;
            } else if (targetNode == sourceNode.getParent()) {
                return targetNode;
            }
        }

        /* Finally, the most general case. We go through all incident shapes and keep track of the highest
         * common ancestor we have found so far. For each new node we process, we distinguish three cases:
         *   1. It is the common ancestor itself or one of the common ancestor's descendants. In this case,
         *      nothing needs to be done.
         *   2. It is a sibling of the common ancestor. In this case, the common ancestor is the new node's
         *      parent.
         *   3. The common ancestor is either a descendant of the new shape or the two are not related in
         *      any way. In this case, we simply traverse down the ancestor hierarchy of both nodes and
         *      set the common ancestor to the lowest common ancestor of the two that we find.
         */

        Iterator<ElkConnectableShape> incidentShapes = allIncidentShapes(edge).iterator();
        ElkNode commonAncestor = connectableShapeToNode(incidentShapes.next());

        while (incidentShapes.hasNext()) {
            ElkNode incidentNode = connectableShapeToNode(incidentShapes.next());

            // Check if the current common ancestor is not an ancestor to the new node, in which case we need to act
            if (incidentNode != commonAncestor && !isDescendant(incidentNode, commonAncestor)) {
                if (incidentNode.getParent() == commonAncestor.getParent()) {
                    // The two nodes are siblings, the common ancestor is their parent
                    commonAncestor = incidentNode.getParent();
                } else {
                    commonAncestor = findLowestCommonAncestor(commonAncestor, incidentNode);

                    if (commonAncestor == null) {
                        // The nodes are not part of the same graph. Abort.
                        return null;
                    }
                }
            }
        }

        return commonAncestor;
    }

    /**
     * Returns the lowest common ancestor of the given two nodes. If the two nodes are not part of the same graph
     * (that is, their root nodes differ), there is no common ancestor.
     *
     * @param node1 the first node.
     * @param node2 the second node.
     * @return the lowest common ancestor or {@code null} if there is none.
     */
    public static ElkNode findLowestCommonAncestor(final ElkNode node1, final ElkNode node2) {
        // Retrieve iterators over the node ancestors
        List<ElkNode> ancestors1 = Lists.newArrayList(new AncestorIterator(node1, true));
        ListIterator<ElkNode> iterator1 = ancestors1.listIterator(ancestors1.size());

        List<ElkNode> ancestors2 = Lists.newArrayList(new AncestorIterator(node2, true));
        ListIterator<ElkNode> iterator2 = ancestors2.listIterator(ancestors2.size());

        // Traverse the ancestor hierarchies from the end as longs as the elements we find are the same
        ElkNode commonAncestor = null;

        while (iterator1.hasPrevious() && iterator2.hasPrevious()) {
            ElkNode ancestor1 = iterator1.previous();
            ElkNode ancestor2 = iterator2.previous();

            if (ancestor1 == ancestor2) {
                commonAncestor = ancestor1;
            } else {
                // The ancestral lines differ; no need to continue
                break;
            }
        }

        return commonAncestor;
    }



    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Convenience

    /**
     * Returns an iterable which contains all edges incoming to a node, whether directly connected or
     * through a port.
     *
     * @param node the node whose incoming edges to gather.
     * @return iterable with all incoming edges.
     */
    public static Iterable<ElkEdge> allIncomingEdges(final ElkNode node) {
        List<Iterable<ElkEdge>> incomingEdgeIterables = Lists.newArrayListWithCapacity(1 + node.getPorts().size());

        incomingEdgeIterables.add(node.getIncomingEdges());
        for (ElkPort port : node.getPorts()) {
            incomingEdgeIterables.add(port.getIncomingEdges());
        }

        return Iterables.concat(incomingEdgeIterables);
    }

    /**
     * Returns an iterable which contains all edges outgoing from a node, whether directly connected or
     * through a port.
     *
     * @param node the node whose outgoing edges to gather.
     * @return iterable with all outgoing edges.
     */
    public static Iterable<ElkEdge> allOutgoingEdges(final ElkNode node) {
        List<Iterable<ElkEdge>> outgoingEdgeIterables = Lists.newArrayListWithCapacity(1 + node.getPorts().size());

        outgoingEdgeIterables.add(node.getOutgoingEdges());
        for (ElkPort port : node.getPorts()) {
            outgoingEdgeIterables.add(port.getOutgoingEdges());
        }

        return Iterables.concat(outgoingEdgeIterables);
    }

    /**
     * Returns an iterable which contains all edges the given connectable shape is incident to.
     *
     * @param shape the connectable shape whose incident edges to return.
     * @return iterable with all incident edges.
     */
    public static Iterable<ElkEdge> allIncidentEdges(final ElkConnectableShape shape) {
        return Iterables.concat(shape.getIncomingEdges(), shape.getOutgoingEdges());
    }

    /**
     * Returns an iterable which contains all edges incident to a node, whether directly connected or
     * through a port.
     *
     * @param node the node whose incident edges to gather.
     * @return iterable with all incident edges.
     */
    public static Iterable<ElkEdge> allIncidentEdges(final ElkNode node) {
        return Iterables.concat(allOutgoingEdges(node), allIncomingEdges(node));
    }

    /**
     * Returns an iterable which contains all connectable shapes the given edge is incident to.
     *
     * @param edge the edge whose end points to return.
     * @return iterable with all incident shapes.
     */
    public static Iterable<ElkConnectableShape> allIncidentShapes(final ElkEdge edge) {
        return Iterables.concat(edge.getSources(), edge.getTargets());
    }

    /**
     * Gathers all edge sections incident to a given edge section, regardless of whether the sections are
     * incoming or outgoing.
     *
     * @param section the section whose incident sections to gather.
     * @return iterable over all incident sections.
     */
    public static Iterable<ElkEdgeSection> allIncidentSections(final ElkEdgeSection section) {
        return Iterables.concat(section.getIncomingSections(), section.getOutgoingSections());
    }

    /**
     * Determines whether the given child node is a descendant of the given ancestor. This method is not reflexive (a
     * node is not its own descendant).
     *
     * @param child a child node.
     * @param ancestor a prospective ancestory node.
     * @return {@code true} if {@code child} is a direct or indirect child of {@code ancestor}.
     */
    public static boolean isDescendant(final ElkNode child, final ElkNode ancestor) {
        // Go up the hierarchy and see if we find the ancestor
        ElkNode current = child;
        while (current.getParent() != null) {
            current = current.getParent();
            if (current == ancestor) {
                return true;
            }
        }

        // Reached the root node without finding the ancestor
        return false;
    }

    /**
     * Returns the node that represents the graph the given element is part of. The method can return {@code null}.
     * This either indicates that the given element is itself the root node of a graph, or that the element is part
     * of an invalid graph structure.
     *
     * @param element the element whose graph to return.
     * @return the graph that contains the element, or {@code null} if such a graph could not be found.
     */
    public static ElkNode containingGraph(final ElkGraphElement element) {
        ElkGraphElement current = element;

        while (current != null) {
            if (current instanceof ElkEdge) {
                // Edges have a direct reference to their containing graph
                return ((ElkEdge) current).getContainingNode();
            } else if (current instanceof ElkNode) {
                // A node's container is its graph
                return ((ElkNode) current).getParent();
            } else if (current.eContainer() instanceof ElkGraphElement) {
                // Walk up the hierarchy
                current = (ElkGraphElement) current.eContainer();
            } else {
                // Something's wrong
                return null;
            }
        }

        // --> current is null
        return null;
    }

    /**
     * Returns the node that belongs to the given connectable shape. That is, if the shape is a node, that itself is
     * returned. If it is a port, the port's parent node is returned. This method may well return {@code null} if the
     * connectable shape is a port that does not belong to a node.
     *
     * @param connectableShape the shape whose node to return.
     * @return the node that belongs to the shape.
     * @throws NullPointerException if {@code connectableShape} is {@code null}.
     */
    public static ElkNode connectableShapeToNode(final ElkConnectableShape connectableShape) {
        if (connectableShape instanceof ElkNode) {
            return (ElkNode) connectableShape;
        } else if (connectableShape instanceof ElkPort) {
            return ((ElkPort) connectableShape).getParent();
        } else if (connectableShape == null) {
            throw new NullPointerException("connectableShape cannot be null");
        } else {
            // In case the meta model is changed in the distant future...
            throw new UnsupportedOperationException("Only support nodes and ports.");
        }
    }

    /**
     * Returns the port that belongs to the given connectable shape, if any. That is, if the shape is a port, that
     * itself is returned. If it is not, {@code null} is returned.
     *
     * @param connectableShape the shape whose port to return.
     * @return the port that belongs to the shape or {@code null}.
     * @throws NullPointerException if {@code connectableShape} is {@code null}.
     */
    public static ElkPort connectableShapeToPort(final ElkConnectableShape connectableShape) {
        if (connectableShape instanceof ElkPort) {
            return (ElkPort) connectableShape;
        } else if (connectableShape == null) {
            throw new NullPointerException("connectableShape cannot be null");
        } else {
            return null;
        }
    }

    /**
     * Returns the source node of the passed simple edge. That is, the passed edge must have
     * exactly one source and one target.
     *
     * @param simpleEdge the edge whose source node to return.
     * @return the source node of the edge, regardless whether the actual source is a port or a node.
     * @throws IllegalArgumentException if {@code simpleEdge} is not 'simple'.
     */
    public static ElkNode getSourceNode(final ElkEdge simpleEdge) {
        if (simpleEdge.getSources().size() != 1 || simpleEdge.getTargets().size() != 1) {
            throw new IllegalArgumentException("Passed edge is not 'simple'.");
        }
        return connectableShapeToNode(simpleEdge.getSources().get(0));
    }

    /**
     * Returns the target node of the passed simple edge. That is, the passed edge must have
     * exactly one source and one target.
     *
     * @param simpleEdge the edge whose target node to return.
     * @return the source node of the edge, regardless whether the actual source is a port or a node.
     * @throws IllegalArgumentException if {@code simpleEdge} is not 'simple'.
     */
    public static ElkNode getTargetNode(final ElkEdge simpleEdge) {
        if (simpleEdge.getSources().size() != 1 || simpleEdge.getTargets().size() != 1) {
            throw new IllegalArgumentException("Passed edge is not 'simple'.");
        }
        return connectableShapeToNode(simpleEdge.getTargets().get(0));
    }

    /**
     * Returns the source port of the passed simple edge. That is, the passed edge must have
     * exactly one source and one target.
     *
     * @param simpleEdge the edge whose source port to return.
     * @return the source port of the edge. If the edge connects directly to a node, {@code null} is returned.
     * @throws IllegalArgumentException if {@code simpleEdge} is not 'simple'.
     */
    public static ElkPort getSourcePort(final ElkEdge simpleEdge) {
        if (simpleEdge.getSources().size() != 1 || simpleEdge.getTargets().size() != 1) {
            throw new IllegalArgumentException("Passed edge is not 'simple'.");
        }
        return connectableShapeToPort(simpleEdge.getSources().get(0));
    }

    /**
     * Returns the target port of the passed simple edge. That is, the passed edge must have
     * exactly one source and one target.
     *
     * @param simpleEdge the edge whose target port to return.
     * @return the source node of the edge. If the edge connects directly to a node, {@code null} is returned.
     * @throws IllegalArgumentException if {@code simpleEdge} is not 'simple'.
     */
    public static ElkPort getTargetPort(final ElkEdge simpleEdge) {
        if (simpleEdge.getSources().size() != 1 || simpleEdge.getTargets().size() != 1) {
            throw new IllegalArgumentException("Passed edge is not 'simple'.");
        }
        return connectableShapeToPort(simpleEdge.getTargets().get(0));
    }
    
    /**
     * Finds the non-label element that the given label labels. Since labels can be nested, we traverse up the
     * containment hierarchy until we find a non-label element.
     * 
     * @param label
     *            the label whose labeled element to find.
     * @return first non-label ancestor or {@code null} if there is none.
     */
    public static ElkGraphElement elementLabeledBy(final ElkLabel label) {
        EObject element = label.getParent();
        while (element instanceof ElkLabel) {
            element = element.eContainer();
        }
        
        // Ensure that the element we've found is actually an ElkGraphElement
        return element instanceof ElkGraphElement
                ? (ElkGraphElement) element
                : null;
    }
    
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Iteration

    /**
     * A tree iterator that skips properties of {@link EMapPropertyHolder}s.
     */
    private static class PropertiesSkippingTreeIterator extends AbstractTreeIterator<EObject> {
        /** Bogus serial version ID. */
        private static final long serialVersionUID = 1L;


        PropertiesSkippingTreeIterator(final Object object, final boolean includeRoot) {
            super(object, includeRoot);
        }


        @Override
        protected Iterator<EObject> getChildren(final Object object) {
            // We know that the object is an EObject; get an iterator over its content
            Iterator<EObject> iterator = ((EObject) object).eContents().iterator();

            // The iterator will usually be a FeatureIteratorImpl that we can set a feature filter on
            if (iterator instanceof FeatureIteratorImpl) {
                ((FeatureIteratorImpl<EObject>) iterator).filter(new FeatureFilter() {
                    public boolean isIncluded(final EStructuralFeature eStructuralFeature) {
                        // We include everything but properties (layout options)
                        if (eStructuralFeature.getContainerClass().equals(EMapPropertyHolder.class)) {
                            return eStructuralFeature.getFeatureID()
                                    != ElkGraphPackage.EMAP_PROPERTY_HOLDER__PROPERTIES;
                        } else {
                            return true;
                        }
                    }
                });
            }

            return iterator;
        }
    }

    /**
     * Returns an iterator over the EMF tree rooted at the given object which skips properties. This is usually useful
     * for iterating over an ELK graph while modifying the properties of elements, which would otherwise be prone to
     * throwing {@link ConcurrentModificationException}s.
     *
     * @param root the EMF tree's root.
     * @param includeRoot {@code true} if the first returned element should be the tree's root itself.
     * @return the requested iterator.
     */
    public static Iterator<EObject> propertiesSkippingIteratorFor(final EObject root, final boolean includeRoot) {
        return new PropertiesSkippingTreeIterator(root, includeRoot);
    }


    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Property values

    /**
     * Returns whether an enumeration value, usually the value of layout option, is considered <em>advanced</em>.
     * The user should have extended knowledge about the effects and underlying methods that are
     * activated when using the layout option value.
     *
     * @param enumValue an enumeration value, usually representing the value of a layout option.
     * @return whether the value is annotated as {@link AdvancedPropertyValue}.
     * @see AdvancedPropertyValue
     */
    public static boolean isAdvancedPropertyValue(final Enum<?> enumValue) {
        // elkjs-exclude-start
        if (enumValue != null) {
            try {
                Annotation[] annotations =
                        enumValue.getClass().getField(enumValue.name()).getAnnotations();
                return Arrays.stream(annotations)
                        .anyMatch(a -> a instanceof AdvancedPropertyValue);
            } catch (NoSuchFieldException | SecurityException e) {
                return false;
            }
        }
        // elkjs-exclude-end

        return false;
    }

    /**
     * Returns whether an enumeration value, usually the value of layout option, is considered <em>experimental</em>.
     * Experimental layout options may result in unexpected behavior, for instance because their implementation has not
     * been thoroughly tested yet.
     *
     * @param enumValue
     *            an enumeration value, usually representing the value of a layout option.
     * @return whether the value is annotated as {@link ExperimentalPropertyValue}.
     * @see ExperimentalPropertyValue
     */
    public static boolean isExperimentalPropertyValue(final Enum<?> enumValue) {
        // elkjs-exclude-start
        if (enumValue != null) {
            try {
                Annotation[] annotations =
                        enumValue.getClass().getField(enumValue.name()).getAnnotations();
                return Arrays.stream(annotations)
                        .anyMatch(a -> a instanceof ExperimentalPropertyValue);
            } catch (NoSuchFieldException | SecurityException e) {
                return false;
            }
        }
        // elkjs-exclude-end

        return false;
    }

    

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Privates

    /**
     * Private constructor, don't call.
     */
    private ElkGraphUtil() {
        // Nothing to do here...
    }


    /**
     * An iterator which, given a starting node, walks up the ancestor tree. The iteration can either start with the
     * node itself or with its parent.
     */
    private static class AncestorIterator implements Iterator<ElkNode> {

        /** The next node we will return. */
        private ElkNode nextNode;


        /**
         * Creates a new iterator.
         *
         * @param startNode the node whose ancestors we want to travel along.
         * @param includeNode {@code true} if {@code startNode} should be the first thing we return, {@code false}
         *                    if its parent should be the first thing.
         */
        AncestorIterator(final ElkNode startNode, final boolean includeNode) {
            nextNode = includeNode ? startNode : startNode.getParent();
        }


        /* (non-Javadoc)
         * @see java.util.Iterator#hasNext()
         */
        @Override
        public boolean hasNext() {
            return nextNode != null;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#next()
         */
        @Override
        public ElkNode next() {
            if (nextNode == null) {
                throw new NoSuchElementException("There is no more element.");
            }
            ElkNode next = nextNode;
            nextNode = nextNode.getParent();
            return next;
        }

    }

}
