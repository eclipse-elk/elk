/**
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 */
package org.eclipse.elk.graph.impl;

import java.util.Collection;

import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkGraphPackage;

import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.common.util.EList;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.util.EObjectWithInverseResolvingEList;
import org.eclipse.emf.ecore.util.InternalEList;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Elk Connectable Shape</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkConnectableShapeImpl#getOutgoingEdges <em>Outgoing Edges</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkConnectableShapeImpl#getIncomingEdges <em>Incoming Edges</em>}</li>
 * </ul>
 *
 * @generated
 */
public abstract class ElkConnectableShapeImpl extends ElkShapeImpl implements ElkConnectableShape {
    /**
     * The cached value of the '{@link #getOutgoingEdges() <em>Outgoing Edges</em>}' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getOutgoingEdges()
     * @generated
     * @ordered
     */
    protected EList<ElkEdge> outgoingEdges;

    /**
     * The cached value of the '{@link #getIncomingEdges() <em>Incoming Edges</em>}' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getIncomingEdges()
     * @generated
     * @ordered
     */
    protected EList<ElkEdge> incomingEdges;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    protected ElkConnectableShapeImpl() {
        super();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    protected EClass eStaticClass() {
        return ElkGraphPackage.Literals.ELK_CONNECTABLE_SHAPE;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkEdge> getOutgoingEdges() {
        if (outgoingEdges == null) {
            outgoingEdges = new EObjectWithInverseResolvingEList.ManyInverse<ElkEdge>(ElkEdge.class, this, ElkGraphPackage.ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES, ElkGraphPackage.ELK_EDGE__SOURCES);
        }
        return outgoingEdges;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkEdge> getIncomingEdges() {
        if (incomingEdges == null) {
            incomingEdges = new EObjectWithInverseResolvingEList.ManyInverse<ElkEdge>(ElkEdge.class, this, ElkGraphPackage.ELK_CONNECTABLE_SHAPE__INCOMING_EDGES, ElkGraphPackage.ELK_EDGE__TARGETS);
        }
        return incomingEdges;
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
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES:
                return ((InternalEList<InternalEObject>)(InternalEList<?>)getOutgoingEdges()).basicAdd(otherEnd, msgs);
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__INCOMING_EDGES:
                return ((InternalEList<InternalEObject>)(InternalEList<?>)getIncomingEdges()).basicAdd(otherEnd, msgs);
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
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES:
                return ((InternalEList<?>)getOutgoingEdges()).basicRemove(otherEnd, msgs);
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__INCOMING_EDGES:
                return ((InternalEList<?>)getIncomingEdges()).basicRemove(otherEnd, msgs);
        }
        return super.eInverseRemove(otherEnd, featureID, msgs);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public Object eGet(int featureID, boolean resolve, boolean coreType) {
        switch (featureID) {
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES:
                return getOutgoingEdges();
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__INCOMING_EDGES:
                return getIncomingEdges();
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
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES:
                getOutgoingEdges().clear();
                getOutgoingEdges().addAll((Collection<? extends ElkEdge>)newValue);
                return;
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__INCOMING_EDGES:
                getIncomingEdges().clear();
                getIncomingEdges().addAll((Collection<? extends ElkEdge>)newValue);
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
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES:
                getOutgoingEdges().clear();
                return;
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__INCOMING_EDGES:
                getIncomingEdges().clear();
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
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES:
                return outgoingEdges != null && !outgoingEdges.isEmpty();
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE__INCOMING_EDGES:
                return incomingEdges != null && !incomingEdges.isEmpty();
        }
        return super.eIsSet(featureID);
    }

} //ElkConnectableShapeImpl
