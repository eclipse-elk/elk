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

import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkGraphPackage;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.common.util.EList;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;

import org.eclipse.emf.ecore.util.EObjectContainmentWithInverseEList;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.emf.ecore.util.InternalEList;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Elk Node</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkNodeImpl#getPorts <em>Ports</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkNodeImpl#getChildren <em>Children</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkNodeImpl#getParent <em>Parent</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkNodeImpl#getContainedEdges <em>Contained Edges</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkNodeImpl#isHierarchical <em>Hierarchical</em>}</li>
 * </ul>
 *
 * @generated
 */
public class ElkNodeImpl extends ElkConnectableShapeImpl implements ElkNode {
    /**
     * The cached value of the '{@link #getPorts() <em>Ports</em>}' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getPorts()
     * @generated
     * @ordered
     */
    protected EList<ElkPort> ports;

    /**
     * The cached value of the '{@link #getChildren() <em>Children</em>}' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getChildren()
     * @generated
     * @ordered
     */
    protected EList<ElkNode> children;

    /**
     * The cached value of the '{@link #getContainedEdges() <em>Contained Edges</em>}' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getContainedEdges()
     * @generated
     * @ordered
     */
    protected EList<ElkEdge> containedEdges;

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
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    protected ElkNodeImpl() {
        super();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    protected EClass eStaticClass() {
        return ElkGraphPackage.Literals.ELK_NODE;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkPort> getPorts() {
        if (ports == null) {
            ports = new EObjectContainmentWithInverseEList<ElkPort>(ElkPort.class, this, ElkGraphPackage.ELK_NODE__PORTS, ElkGraphPackage.ELK_PORT__PARENT);
        }
        return ports;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkNode> getChildren() {
        if (children == null) {
            children = new EObjectContainmentWithInverseEList<ElkNode>(ElkNode.class, this, ElkGraphPackage.ELK_NODE__CHILDREN, ElkGraphPackage.ELK_NODE__PARENT);
        }
        return children;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkNode getParent() {
        if (eContainerFeatureID() != ElkGraphPackage.ELK_NODE__PARENT) return null;
        return (ElkNode)eInternalContainer();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public NotificationChain basicSetParent(ElkNode newParent, NotificationChain msgs) {
        msgs = eBasicSetContainer((InternalEObject)newParent, ElkGraphPackage.ELK_NODE__PARENT, msgs);
        return msgs;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setParent(ElkNode newParent) {
        if (newParent != eInternalContainer() || (eContainerFeatureID() != ElkGraphPackage.ELK_NODE__PARENT && newParent != null)) {
            if (EcoreUtil.isAncestor(this, newParent))
                throw new IllegalArgumentException("Recursive containment not allowed for " + toString());
            NotificationChain msgs = null;
            if (eInternalContainer() != null)
                msgs = eBasicRemoveFromContainer(msgs);
            if (newParent != null)
                msgs = ((InternalEObject)newParent).eInverseAdd(this, ElkGraphPackage.ELK_NODE__CHILDREN, ElkNode.class, msgs);
            msgs = basicSetParent(newParent, msgs);
            if (msgs != null) msgs.dispatch();
        }
        else if (eNotificationRequired())
            eNotify(new ENotificationImpl(this, Notification.SET, ElkGraphPackage.ELK_NODE__PARENT, newParent, newParent));
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkEdge> getContainedEdges() {
        if (containedEdges == null) {
            containedEdges = new EObjectContainmentWithInverseEList<ElkEdge>(ElkEdge.class, this, ElkGraphPackage.ELK_NODE__CONTAINED_EDGES, ElkGraphPackage.ELK_EDGE__CONTAINING_NODE);
        }
        return containedEdges;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public boolean isHierarchical() {
        return getChildren().size() > 0;
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
            case ElkGraphPackage.ELK_NODE__PORTS:
                return ((InternalEList<InternalEObject>)(InternalEList<?>)getPorts()).basicAdd(otherEnd, msgs);
            case ElkGraphPackage.ELK_NODE__CHILDREN:
                return ((InternalEList<InternalEObject>)(InternalEList<?>)getChildren()).basicAdd(otherEnd, msgs);
            case ElkGraphPackage.ELK_NODE__PARENT:
                if (eInternalContainer() != null)
                    msgs = eBasicRemoveFromContainer(msgs);
                return basicSetParent((ElkNode)otherEnd, msgs);
            case ElkGraphPackage.ELK_NODE__CONTAINED_EDGES:
                return ((InternalEList<InternalEObject>)(InternalEList<?>)getContainedEdges()).basicAdd(otherEnd, msgs);
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
            case ElkGraphPackage.ELK_NODE__PORTS:
                return ((InternalEList<?>)getPorts()).basicRemove(otherEnd, msgs);
            case ElkGraphPackage.ELK_NODE__CHILDREN:
                return ((InternalEList<?>)getChildren()).basicRemove(otherEnd, msgs);
            case ElkGraphPackage.ELK_NODE__PARENT:
                return basicSetParent(null, msgs);
            case ElkGraphPackage.ELK_NODE__CONTAINED_EDGES:
                return ((InternalEList<?>)getContainedEdges()).basicRemove(otherEnd, msgs);
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
            case ElkGraphPackage.ELK_NODE__PARENT:
                return eInternalContainer().eInverseRemove(this, ElkGraphPackage.ELK_NODE__CHILDREN, ElkNode.class, msgs);
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
            case ElkGraphPackage.ELK_NODE__PORTS:
                return getPorts();
            case ElkGraphPackage.ELK_NODE__CHILDREN:
                return getChildren();
            case ElkGraphPackage.ELK_NODE__PARENT:
                return getParent();
            case ElkGraphPackage.ELK_NODE__CONTAINED_EDGES:
                return getContainedEdges();
            case ElkGraphPackage.ELK_NODE__HIERARCHICAL:
                return isHierarchical();
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
            case ElkGraphPackage.ELK_NODE__PORTS:
                getPorts().clear();
                getPorts().addAll((Collection<? extends ElkPort>)newValue);
                return;
            case ElkGraphPackage.ELK_NODE__CHILDREN:
                getChildren().clear();
                getChildren().addAll((Collection<? extends ElkNode>)newValue);
                return;
            case ElkGraphPackage.ELK_NODE__PARENT:
                setParent((ElkNode)newValue);
                return;
            case ElkGraphPackage.ELK_NODE__CONTAINED_EDGES:
                getContainedEdges().clear();
                getContainedEdges().addAll((Collection<? extends ElkEdge>)newValue);
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
            case ElkGraphPackage.ELK_NODE__PORTS:
                getPorts().clear();
                return;
            case ElkGraphPackage.ELK_NODE__CHILDREN:
                getChildren().clear();
                return;
            case ElkGraphPackage.ELK_NODE__PARENT:
                setParent((ElkNode)null);
                return;
            case ElkGraphPackage.ELK_NODE__CONTAINED_EDGES:
                getContainedEdges().clear();
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
            case ElkGraphPackage.ELK_NODE__PORTS:
                return ports != null && !ports.isEmpty();
            case ElkGraphPackage.ELK_NODE__CHILDREN:
                return children != null && !children.isEmpty();
            case ElkGraphPackage.ELK_NODE__PARENT:
                return getParent() != null;
            case ElkGraphPackage.ELK_NODE__CONTAINED_EDGES:
                return containedEdges != null && !containedEdges.isEmpty();
            case ElkGraphPackage.ELK_NODE__HIERARCHICAL:
                return isHierarchical() != HIERARCHICAL_EDEFAULT;
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
        
        StringBuilder builder = new StringBuilder("ElkNode");
        // Use identifier or labels
        String id = getIdentifier();
        if (id != null && !id.isEmpty()) {
            builder.append(" \"").append(id).append("\"");
        } else if (getLabels().size() > 0) {
            String text = getLabels().get(0).getText();
            if (text != null && !text.isEmpty()) {
                builder.append(" \"").append(text).append("\"");
            }
        }
        // position and dimension
        builder.append(" (").append(x).append(",").append(y).append(" | ")
               .append(width).append(",").append(height).append(")");

        return builder.toString();
    }

} //ElkNodeImpl
