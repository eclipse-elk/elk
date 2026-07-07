/**
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 */
package org.eclipse.elk.graph;

import org.eclipse.emf.common.util.EList;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Elk Edge Section</b></em>'.
 * <!-- end-user-doc -->
 *
 * <!-- begin-model-doc -->
 * An edge section defines the routing of an edge, or a part of that routing. If the edge is a regular edge (as opposed to a hyperedge), it will have a single edge section that connects to the edge's source element and target element. The section will then completely define the edge's start point, its end point, and all of its bend points. This is a special case of the more general case described below.
 * 
 * <p>If the section's parent edge is a hyperedge, defining the routing will be more complicated. There will be enough edge sections to connect all of the edge's souce and target elements. The sections will effectively define a <em>routing graph</em>: all sections in the graph will connect to other sections and/or sources and targets of the edge, each effectively defining a part of the complex route the edge will take. We call an edge section that connects to at least one {@link ElkConnectableShape} an <em>outer section</em>. Edge sections that connect only to other edge sections are referred to as <em>inner sections</em>.</p>
 * 
 * <p>Conceptually, the routing graph would be undirected. The way references work in EMF, however, forces us to distinguish between a section's incoming and outgoing sections. This is not much of a problem, however: each routing graph can be made acyclic.</p>
 * 
 * <p>All coordinates that define a section's route are relative to the origin of its edge's containing node.</p>
 * 
 * <p>Note that edge sections are property holders to allow algorithms to pass more detailed information about an edge section back to the client.</p>
 * <!-- end-model-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getStartX <em>Start X</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getStartY <em>Start Y</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getEndX <em>End X</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getEndY <em>End Y</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getBendPoints <em>Bend Points</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getParent <em>Parent</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getOutgoingShape <em>Outgoing Shape</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getIncomingShape <em>Incoming Shape</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getOutgoingSections <em>Outgoing Sections</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getIncomingSections <em>Incoming Sections</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.ElkEdgeSection#getIdentifier <em>Identifier</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection()
 * @model
 * @generated
 */
public interface ElkEdgeSection extends EMapPropertyHolder {
    /**
     * Returns the value of the '<em><b>Start X</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * X coordinate of the section's start point, relative to the origin of the edge's containing node.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Start X</em>' attribute.
     * @see #setStartX(double)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_StartX()
     * @model
     * @generated
     */
    double getStartX();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkEdgeSection#getStartX <em>Start X</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Start X</em>' attribute.
     * @see #getStartX()
     * @generated
     */
    void setStartX(double value);

    /**
     * Returns the value of the '<em><b>Start Y</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Y coordinate of the section's start point, relative to the origin of the edge's containing node.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Start Y</em>' attribute.
     * @see #setStartY(double)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_StartY()
     * @model
     * @generated
     */
    double getStartY();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkEdgeSection#getStartY <em>Start Y</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Start Y</em>' attribute.
     * @see #getStartY()
     * @generated
     */
    void setStartY(double value);

    /**
     * Returns the value of the '<em><b>End X</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * X coordinate of the section's end point, relative to the origin of the edge's containing node.
     * <!-- end-model-doc -->
     * @return the value of the '<em>End X</em>' attribute.
     * @see #setEndX(double)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_EndX()
     * @model
     * @generated
     */
    double getEndX();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkEdgeSection#getEndX <em>End X</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>End X</em>' attribute.
     * @see #getEndX()
     * @generated
     */
    void setEndX(double value);

    /**
     * Returns the value of the '<em><b>End Y</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Y coordinate of the section's end point, relative to the origin of the edge's containing node.
     * <!-- end-model-doc -->
     * @return the value of the '<em>End Y</em>' attribute.
     * @see #setEndY(double)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_EndY()
     * @model
     * @generated
     */
    double getEndY();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkEdgeSection#getEndY <em>End Y</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>End Y</em>' attribute.
     * @see #getEndY()
     * @generated
     */
    void setEndY(double value);

    /**
     * Returns the value of the '<em><b>Bend Points</b></em>' containment reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkBendPoint}.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The section's list of bend points. May well be empty if the section represents a straight line.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Bend Points</em>' containment reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_BendPoints()
     * @model containment="true"
     * @generated
     */
    EList<ElkBendPoint> getBendPoints();

    /**
     * Returns the value of the '<em><b>Parent</b></em>' container reference.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkEdge#getSections <em>Sections</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The edge this section belongs to.
     * 
     * <p>Setting the parent edge automatically updates its list of edge sections.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Parent</em>' container reference.
     * @see #setParent(ElkEdge)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_Parent()
     * @see org.eclipse.elk.graph.ElkEdge#getSections
     * @model opposite="sections" transient="false"
     * @generated
     */
    ElkEdge getParent();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkEdgeSection#getParent <em>Parent</em>}' container reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Parent</em>' container reference.
     * @see #getParent()
     * @generated
     */
    void setParent(ElkEdge value);

    /**
     * Returns the value of the '<em><b>Outgoing Shape</b></em>' reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The shape this section ends at, if any. If there is one, this section is an <em>outer section</em>.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Outgoing Shape</em>' reference.
     * @see #setOutgoingShape(ElkConnectableShape)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_OutgoingShape()
     * @model
     * @generated
     */
    ElkConnectableShape getOutgoingShape();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkEdgeSection#getOutgoingShape <em>Outgoing Shape</em>}' reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Outgoing Shape</em>' reference.
     * @see #getOutgoingShape()
     * @generated
     */
    void setOutgoingShape(ElkConnectableShape value);

    /**
     * Returns the value of the '<em><b>Incoming Shape</b></em>' reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * The shape this section starts at, if any. If there is one, this section is an <em>outer section</em>.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Incoming Shape</em>' reference.
     * @see #setIncomingShape(ElkConnectableShape)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_IncomingShape()
     * @model
     * @generated
     */
    ElkConnectableShape getIncomingShape();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkEdgeSection#getIncomingShape <em>Incoming Shape</em>}' reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Incoming Shape</em>' reference.
     * @see #getIncomingShape()
     * @generated
     */
    void setIncomingShape(ElkConnectableShape value);

    /**
     * Returns the value of the '<em><b>Outgoing Sections</b></em>' reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkEdgeSection}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkEdgeSection#getIncomingSections <em>Incoming Sections</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * List of outgoing sections this section is connected to. Must not be empty if this section is an <em>inner section</em> (not connected to a shape).
     * 
     * <p>Adding or removing a section to/from this list automatically updates its list of incoming sections.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Outgoing Sections</em>' reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_OutgoingSections()
     * @see org.eclipse.elk.graph.ElkEdgeSection#getIncomingSections
     * @model opposite="incomingSections"
     * @generated
     */
    EList<ElkEdgeSection> getOutgoingSections();

    /**
     * Returns the value of the '<em><b>Incoming Sections</b></em>' reference list.
     * The list contents are of type {@link org.eclipse.elk.graph.ElkEdgeSection}.
     * It is bidirectional and its opposite is '{@link org.eclipse.elk.graph.ElkEdgeSection#getOutgoingSections <em>Outgoing Sections</em>}'.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * List of incoming sections this section is connected to. Must not be empty if this section is an <em>inner section</em> (not connected to a shape).
     * 
     * <p>Adding or removing a section to/from this list automatically updates its list of outgoing sections.</p>
     * <!-- end-model-doc -->
     * @return the value of the '<em>Incoming Sections</em>' reference list.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_IncomingSections()
     * @see org.eclipse.elk.graph.ElkEdgeSection#getOutgoingSections
     * @model opposite="outgoingSections"
     * @generated
     */
    EList<ElkEdgeSection> getIncomingSections();

    /**
     * Returns the value of the '<em><b>Identifier</b></em>' attribute.
     * <!-- begin-user-doc -->
     * <p>
     * If the meaning of the '<em>Identifier</em>' attribute isn't clear,
     * there really should be more of a description here...
     * </p>
     * <!-- end-user-doc -->
     * @return the value of the '<em>Identifier</em>' attribute.
     * @see #setIdentifier(String)
     * @see org.eclipse.elk.graph.ElkGraphPackage#getElkEdgeSection_Identifier()
     * @model
     * @generated
     */
    String getIdentifier();

    /**
     * Sets the value of the '{@link org.eclipse.elk.graph.ElkEdgeSection#getIdentifier <em>Identifier</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param value the new value of the '<em>Identifier</em>' attribute.
     * @see #getIdentifier()
     * @generated
     */
    void setIdentifier(String value);

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Sets the x and y coordinates of the section's start point simultaneously by calling their respective set methods.
     * <!-- end-model-doc -->
     * @model annotation="http://www.eclipse.org/emf/2002/GenModel body='setStartX(x);\nsetStartY(y);'"
     * @generated
     */
    void setStartLocation(double x, double y);

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Sets the x and y coordinates of the section's end point simultaneously by calling their respective set methods.
     * <!-- end-model-doc -->
     * @model annotation="http://www.eclipse.org/emf/2002/GenModel body='setEndX(x);\nsetEndY(y);'"
     * @generated
     */
    void setEndLocation(double x, double y);

} // ElkEdgeSection
