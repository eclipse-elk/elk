/**
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 */
package org.eclipse.elk.graph.impl;

import java.util.Collection;

import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphPackage;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;
import org.eclipse.emf.common.util.EList;
import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;
import org.eclipse.emf.ecore.impl.ENotificationImpl;
import org.eclipse.emf.ecore.util.EObjectContainmentWithInverseEList;
import org.eclipse.emf.ecore.util.EObjectWithInverseResolvingEList;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.emf.ecore.util.InternalEList;

import com.google.common.base.Joiner;
import com.google.common.base.Strings;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Elk Edge</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeImpl#getContainingNode <em>Containing Node</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeImpl#getSources <em>Sources</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeImpl#getTargets <em>Targets</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeImpl#getSections <em>Sections</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeImpl#isHyperedge <em>Hyperedge</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeImpl#isHierarchical <em>Hierarchical</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeImpl#isSelfloop <em>Selfloop</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeImpl#isConnected <em>Connected</em>}</li>
 * </ul>
 *
 * @generated
 */
public class ElkEdgeImpl extends ElkGraphElementImpl implements ElkEdge {
    /**
     * The cached value of the '{@link #getSources() <em>Sources</em>}' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getSources()
     * @generated
     * @ordered
     */
    protected EList<ElkConnectableShape> sources;

    /**
     * The cached value of the '{@link #getTargets() <em>Targets</em>}' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getTargets()
     * @generated
     * @ordered
     */
    protected EList<ElkConnectableShape> targets;

    /**
     * The cached value of the '{@link #getSections() <em>Sections</em>}' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getSections()
     * @generated
     * @ordered
     */
    protected EList<ElkEdgeSection> sections;

    /**
     * The default value of the '{@link #isHyperedge() <em>Hyperedge</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #isHyperedge()
     * @generated
     * @ordered
     */
    protected static final boolean HYPEREDGE_EDEFAULT = false;

    /**
     * The default value of the '{@link #isHierarchical() <em>Hierarchical</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #isHierarchical()
     * @generated
     * @ordered
     */
    protected static final boolean HIERARCHICAL_EDEFAULT = false;

    /**
     * The default value of the '{@link #isSelfloop() <em>Selfloop</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #isSelfloop()
     * @generated
     * @ordered
     */
    protected static final boolean SELFLOOP_EDEFAULT = false;

    /**
     * The default value of the '{@link #isConnected() <em>Connected</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #isConnected()
     * @generated
     * @ordered
     */
    protected static final boolean CONNECTED_EDEFAULT = false;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    protected ElkEdgeImpl() {
        super();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    protected EClass eStaticClass() {
        return ElkGraphPackage.Literals.ELK_EDGE;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkNode getContainingNode() {
        if (eContainerFeatureID() != ElkGraphPackage.ELK_EDGE__CONTAINING_NODE) return null;
        return (ElkNode)eInternalContainer();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public NotificationChain basicSetContainingNode(ElkNode newContainingNode, NotificationChain msgs) {
        msgs = eBasicSetContainer((InternalEObject)newContainingNode, ElkGraphPackage.ELK_EDGE__CONTAINING_NODE, msgs);
        return msgs;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setContainingNode(ElkNode newContainingNode) {
        if (newContainingNode != eInternalContainer() || (eContainerFeatureID() != ElkGraphPackage.ELK_EDGE__CONTAINING_NODE && newContainingNode != null)) {
            if (EcoreUtil.isAncestor(this, newContainingNode))
                throw new IllegalArgumentException("Recursive containment not allowed for " + toString());
            NotificationChain msgs = null;
            if (eInternalContainer() != null)
                msgs = eBasicRemoveFromContainer(msgs);
            if (newContainingNode != null)
                msgs = ((InternalEObject)newContainingNode).eInverseAdd(this, ElkGraphPackage.ELK_NODE__CONTAINED_EDGES, ElkNode.class, msgs);
            msgs = basicSetContainingNode(newContainingNode, msgs);
            if (msgs != null) msgs.dispatch();
        }
        else if (eNotificationRequired())
            eNotify(new ENotificationImpl(this, Notification.SET, ElkGraphPackage.ELK_EDGE__CONTAINING_NODE, newContainingNode, newContainingNode));
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkConnectableShape> getSources() {
        if (sources == null) {
            sources = new EObjectWithInverseResolvingEList.ManyInverse<ElkConnectableShape>(ElkConnectableShape.class, this, ElkGraphPackage.ELK_EDGE__SOURCES, ElkGraphPackage.ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES);
        }
        return sources;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkConnectableShape> getTargets() {
        if (targets == null) {
            targets = new EObjectWithInverseResolvingEList.ManyInverse<ElkConnectableShape>(ElkConnectableShape.class, this, ElkGraphPackage.ELK_EDGE__TARGETS, ElkGraphPackage.ELK_CONNECTABLE_SHAPE__INCOMING_EDGES);
        }
        return targets;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkEdgeSection> getSections() {
        if (sections == null) {
            sections = new EObjectContainmentWithInverseEList<ElkEdgeSection>(ElkEdgeSection.class, this, ElkGraphPackage.ELK_EDGE__SECTIONS, ElkGraphPackage.ELK_EDGE_SECTION__PARENT);
        }
        return sections;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public boolean isHyperedge() {
        return !(getSources().size() <= 1 && getTargets().size() <= 1);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public boolean isHierarchical() {
        // We're basically iterating over all sources and targets and check whether they are contained in the same
        // graph. If they are, we're a simple edge. If they are not, we're a hierarchical edge.
        ElkNode commonRepresentingNode = null;
        
        for (ElkConnectableShape incidentShape : ElkGraphUtil.allIncidentShapes(this)) {
            // Compute representing node
            ElkNode shapeNode = ElkGraphUtil.connectableShapeToNode(incidentShape);
            
            if (commonRepresentingNode == null) {
                commonRepresentingNode = shapeNode.getParent();
            } else if (commonRepresentingNode != shapeNode.getParent()) {
                // Different hierarchy levels
                return true;
            }
        }
        
        return false;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public boolean isSelfloop() {
        // We're basically iterating over all sources and targets and check whether they all belong to the same node
        ElkNode commonNode = null;
        
        for (ElkConnectableShape incidentShape : ElkGraphUtil.allIncidentShapes(this)) {
            // Compute representing node
            ElkNode shapeNode = ElkGraphUtil.connectableShapeToNode(incidentShape);
            
            if (commonNode == null) {
                commonNode = shapeNode;
            } else if (commonNode != shapeNode) {
                // Not a self loop
                return false;
            }
        }
        
        return true;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public boolean isConnected() {
        return !getSources().isEmpty() && !getTargets().isEmpty();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @SuppressWarnings("unchecked")
    @Override
    public NotificationChain eInverseAdd(InternalEObject otherEnd, int featureID, NotificationChain msgs) {
        switch (featureID) {
            case ElkGraphPackage.ELK_EDGE__CONTAINING_NODE:
                if (eInternalContainer() != null)
                    msgs = eBasicRemoveFromContainer(msgs);
                return basicSetContainingNode((ElkNode)otherEnd, msgs);
            case ElkGraphPackage.ELK_EDGE__SOURCES:
                return ((InternalEList<InternalEObject>)(InternalEList<?>)getSources()).basicAdd(otherEnd, msgs);
            case ElkGraphPackage.ELK_EDGE__TARGETS:
                return ((InternalEList<InternalEObject>)(InternalEList<?>)getTargets()).basicAdd(otherEnd, msgs);
            case ElkGraphPackage.ELK_EDGE__SECTIONS:
                return ((InternalEList<InternalEObject>)(InternalEList<?>)getSections()).basicAdd(otherEnd, msgs);
        }
        return super.eInverseAdd(otherEnd, featureID, msgs);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public NotificationChain eInverseRemove(InternalEObject otherEnd, int featureID, NotificationChain msgs) {
        switch (featureID) {
            case ElkGraphPackage.ELK_EDGE__CONTAINING_NODE:
                return basicSetContainingNode(null, msgs);
            case ElkGraphPackage.ELK_EDGE__SOURCES:
                return ((InternalEList<?>)getSources()).basicRemove(otherEnd, msgs);
            case ElkGraphPackage.ELK_EDGE__TARGETS:
                return ((InternalEList<?>)getTargets()).basicRemove(otherEnd, msgs);
            case ElkGraphPackage.ELK_EDGE__SECTIONS:
                return ((InternalEList<?>)getSections()).basicRemove(otherEnd, msgs);
        }
        return super.eInverseRemove(otherEnd, featureID, msgs);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public NotificationChain eBasicRemoveFromContainerFeature(NotificationChain msgs) {
        switch (eContainerFeatureID()) {
            case ElkGraphPackage.ELK_EDGE__CONTAINING_NODE:
                return eInternalContainer().eInverseRemove(this, ElkGraphPackage.ELK_NODE__CONTAINED_EDGES, ElkNode.class, msgs);
        }
        return super.eBasicRemoveFromContainerFeature(msgs);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public Object eGet(int featureID, boolean resolve, boolean coreType) {
        switch (featureID) {
            case ElkGraphPackage.ELK_EDGE__CONTAINING_NODE:
                return getContainingNode();
            case ElkGraphPackage.ELK_EDGE__SOURCES:
                return getSources();
            case ElkGraphPackage.ELK_EDGE__TARGETS:
                return getTargets();
            case ElkGraphPackage.ELK_EDGE__SECTIONS:
                return getSections();
            case ElkGraphPackage.ELK_EDGE__HYPEREDGE:
                return isHyperedge();
            case ElkGraphPackage.ELK_EDGE__HIERARCHICAL:
                return isHierarchical();
            case ElkGraphPackage.ELK_EDGE__SELFLOOP:
                return isSelfloop();
            case ElkGraphPackage.ELK_EDGE__CONNECTED:
                return isConnected();
        }
        return super.eGet(featureID, resolve, coreType);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @SuppressWarnings("unchecked")
    @Override
    public void eSet(int featureID, Object newValue) {
        switch (featureID) {
            case ElkGraphPackage.ELK_EDGE__CONTAINING_NODE:
                setContainingNode((ElkNode)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE__SOURCES:
                getSources().clear();
                getSources().addAll((Collection<? extends ElkConnectableShape>)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE__TARGETS:
                getTargets().clear();
                getTargets().addAll((Collection<? extends ElkConnectableShape>)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE__SECTIONS:
                getSections().clear();
                getSections().addAll((Collection<? extends ElkEdgeSection>)newValue);
                return;
        }
        super.eSet(featureID, newValue);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public void eUnset(int featureID) {
        switch (featureID) {
            case ElkGraphPackage.ELK_EDGE__CONTAINING_NODE:
                setContainingNode((ElkNode)null);
                return;
            case ElkGraphPackage.ELK_EDGE__SOURCES:
                getSources().clear();
                return;
            case ElkGraphPackage.ELK_EDGE__TARGETS:
                getTargets().clear();
                return;
            case ElkGraphPackage.ELK_EDGE__SECTIONS:
                getSections().clear();
                return;
        }
        super.eUnset(featureID);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public boolean eIsSet(int featureID) {
        switch (featureID) {
            case ElkGraphPackage.ELK_EDGE__CONTAINING_NODE:
                return getContainingNode() != null;
            case ElkGraphPackage.ELK_EDGE__SOURCES:
                return sources != null && !sources.isEmpty();
            case ElkGraphPackage.ELK_EDGE__TARGETS:
                return targets != null && !targets.isEmpty();
            case ElkGraphPackage.ELK_EDGE__SECTIONS:
                return sections != null && !sections.isEmpty();
            case ElkGraphPackage.ELK_EDGE__HYPEREDGE:
                return isHyperedge() != HYPEREDGE_EDEFAULT;
            case ElkGraphPackage.ELK_EDGE__HIERARCHICAL:
                return isHierarchical() != HIERARCHICAL_EDEFAULT;
            case ElkGraphPackage.ELK_EDGE__SELFLOOP:
                return isSelfloop() != SELFLOOP_EDEFAULT;
            case ElkGraphPackage.ELK_EDGE__CONNECTED:
                return isConnected() != CONNECTED_EDEFAULT;
        }
        return super.eIsSet(featureID);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    @Override
    public String toString() {
        if (eIsProxy()) return super.toString();

        StringBuilder builder = new StringBuilder("ElkEdge");
        // Use identifier or labels
        String id = getIdentifier();
        if (!Strings.isNullOrEmpty(id)) {
            builder.append(" \"").append(id).append("\"");
        } else if (getLabels().size() > 0) {
            String text = getLabels().get(0).getText();
            if (!Strings.isNullOrEmpty(text)) {
                builder.append(" \"").append(text).append("\"");
            }
        }
        
        boolean hyperedge = isHyperedge();
        if (hyperedge) {
            builder.append(" [");
        } else {
            builder.append(" ");
        }
        builder.append(Joiner.on(", ").join(sources));
        if (hyperedge) {
            builder.append("]");
        }
        builder.append(" -> ");
        if (hyperedge) {
            builder.append("[");
        }
        builder.append(Joiner.on(", ").join(targets));
        if (hyperedge) {
            builder.append("]");
        }
        
        return builder.toString();
    }
    
} //ElkEdgeImpl
