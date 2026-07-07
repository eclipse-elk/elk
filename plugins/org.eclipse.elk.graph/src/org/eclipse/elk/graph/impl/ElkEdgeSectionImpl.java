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

import org.eclipse.elk.graph.ElkBendPoint;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphPackage;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.common.util.EList;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;

import org.eclipse.emf.ecore.util.EObjectContainmentEList;
import org.eclipse.emf.ecore.util.EObjectWithInverseResolvingEList;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.emf.ecore.util.InternalEList;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Elk Edge Section</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getStartX <em>Start X</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getStartY <em>Start Y</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getEndX <em>End X</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getEndY <em>End Y</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getBendPoints <em>Bend Points</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getParent <em>Parent</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getOutgoingShape <em>Outgoing Shape</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getIncomingShape <em>Incoming Shape</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getOutgoingSections <em>Outgoing Sections</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getIncomingSections <em>Incoming Sections</em>}</li>
 *   <li>{@link org.eclipse.elk.graph.impl.ElkEdgeSectionImpl#getIdentifier <em>Identifier</em>}</li>
 * </ul>
 *
 * @generated
 */
public class ElkEdgeSectionImpl extends EMapPropertyHolderImpl implements ElkEdgeSection {
    /**
     * The default value of the '{@link #getStartX() <em>Start X</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getStartX()
     * @generated
     * @ordered
     */
    protected static final double START_X_EDEFAULT = 0.0;

    /**
     * The cached value of the '{@link #getStartX() <em>Start X</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getStartX()
     * @generated
     * @ordered
     */
    protected double startX = START_X_EDEFAULT;

    /**
     * The default value of the '{@link #getStartY() <em>Start Y</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getStartY()
     * @generated
     * @ordered
     */
    protected static final double START_Y_EDEFAULT = 0.0;

    /**
     * The cached value of the '{@link #getStartY() <em>Start Y</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getStartY()
     * @generated
     * @ordered
     */
    protected double startY = START_Y_EDEFAULT;

    /**
     * The default value of the '{@link #getEndX() <em>End X</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getEndX()
     * @generated
     * @ordered
     */
    protected static final double END_X_EDEFAULT = 0.0;

    /**
     * The cached value of the '{@link #getEndX() <em>End X</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getEndX()
     * @generated
     * @ordered
     */
    protected double endX = END_X_EDEFAULT;

    /**
     * The default value of the '{@link #getEndY() <em>End Y</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getEndY()
     * @generated
     * @ordered
     */
    protected static final double END_Y_EDEFAULT = 0.0;

    /**
     * The cached value of the '{@link #getEndY() <em>End Y</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getEndY()
     * @generated
     * @ordered
     */
    protected double endY = END_Y_EDEFAULT;

    /**
     * The cached value of the '{@link #getBendPoints() <em>Bend Points</em>}' containment reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getBendPoints()
     * @generated
     * @ordered
     */
    protected EList<ElkBendPoint> bendPoints;

    /**
     * The cached value of the '{@link #getOutgoingShape() <em>Outgoing Shape</em>}' reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getOutgoingShape()
     * @generated
     * @ordered
     */
    protected ElkConnectableShape outgoingShape;

    /**
     * The cached value of the '{@link #getIncomingShape() <em>Incoming Shape</em>}' reference.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getIncomingShape()
     * @generated
     * @ordered
     */
    protected ElkConnectableShape incomingShape;

    /**
     * The cached value of the '{@link #getOutgoingSections() <em>Outgoing Sections</em>}' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getOutgoingSections()
     * @generated
     * @ordered
     */
    protected EList<ElkEdgeSection> outgoingSections;

    /**
     * The cached value of the '{@link #getIncomingSections() <em>Incoming Sections</em>}' reference list.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getIncomingSections()
     * @generated
     * @ordered
     */
    protected EList<ElkEdgeSection> incomingSections;

    /**
     * The default value of the '{@link #getIdentifier() <em>Identifier</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getIdentifier()
     * @generated
     * @ordered
     */
    protected static final String IDENTIFIER_EDEFAULT = null;

    /**
     * The cached value of the '{@link #getIdentifier() <em>Identifier</em>}' attribute.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getIdentifier()
     * @generated
     * @ordered
     */
    protected String identifier = IDENTIFIER_EDEFAULT;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    protected ElkEdgeSectionImpl() {
        super();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    protected EClass eStaticClass() {
        return ElkGraphPackage.Literals.ELK_EDGE_SECTION;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public double getStartX() {
        return startX;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setStartX(double newStartX) {
        double oldStartX = startX;
        startX = newStartX;
        if (eNotificationRequired())
            eNotify(new ENotificationImpl(this, Notification.SET, ElkGraphPackage.ELK_EDGE_SECTION__START_X, oldStartX, startX));
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public double getStartY() {
        return startY;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setStartY(double newStartY) {
        double oldStartY = startY;
        startY = newStartY;
        if (eNotificationRequired())
            eNotify(new ENotificationImpl(this, Notification.SET, ElkGraphPackage.ELK_EDGE_SECTION__START_Y, oldStartY, startY));
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public double getEndX() {
        return endX;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setEndX(double newEndX) {
        double oldEndX = endX;
        endX = newEndX;
        if (eNotificationRequired())
            eNotify(new ENotificationImpl(this, Notification.SET, ElkGraphPackage.ELK_EDGE_SECTION__END_X, oldEndX, endX));
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public double getEndY() {
        return endY;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setEndY(double newEndY) {
        double oldEndY = endY;
        endY = newEndY;
        if (eNotificationRequired())
            eNotify(new ENotificationImpl(this, Notification.SET, ElkGraphPackage.ELK_EDGE_SECTION__END_Y, oldEndY, endY));
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkBendPoint> getBendPoints() {
        if (bendPoints == null) {
            bendPoints = new EObjectContainmentEList<ElkBendPoint>(ElkBendPoint.class, this, ElkGraphPackage.ELK_EDGE_SECTION__BEND_POINTS);
        }
        return bendPoints;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkEdge getParent() {
        if (eContainerFeatureID() != ElkGraphPackage.ELK_EDGE_SECTION__PARENT) return null;
        return (ElkEdge)eInternalContainer();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public NotificationChain basicSetParent(ElkEdge newParent, NotificationChain msgs) {
        msgs = eBasicSetContainer((InternalEObject)newParent, ElkGraphPackage.ELK_EDGE_SECTION__PARENT, msgs);
        return msgs;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setParent(ElkEdge newParent) {
        if (newParent != eInternalContainer() || (eContainerFeatureID() != ElkGraphPackage.ELK_EDGE_SECTION__PARENT && newParent != null)) {
            if (EcoreUtil.isAncestor(this, newParent))
                throw new IllegalArgumentException("Recursive containment not allowed for " + toString());
            NotificationChain msgs = null;
            if (eInternalContainer() != null)
                msgs = eBasicRemoveFromContainer(msgs);
            if (newParent != null)
                msgs = ((InternalEObject)newParent).eInverseAdd(this, ElkGraphPackage.ELK_EDGE__SECTIONS, ElkEdge.class, msgs);
            msgs = basicSetParent(newParent, msgs);
            if (msgs != null) msgs.dispatch();
        }
        else if (eNotificationRequired())
            eNotify(new ENotificationImpl(this, Notification.SET, ElkGraphPackage.ELK_EDGE_SECTION__PARENT, newParent, newParent));
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkConnectableShape getOutgoingShape() {
        if (outgoingShape != null && outgoingShape.eIsProxy()) {
            InternalEObject oldOutgoingShape = (InternalEObject)outgoingShape;
            outgoingShape = (ElkConnectableShape)eResolveProxy(oldOutgoingShape);
            if (outgoingShape != oldOutgoingShape) {
                if (eNotificationRequired())
                    eNotify(new ENotificationImpl(this, Notification.RESOLVE, ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SHAPE, oldOutgoingShape, outgoingShape));
            }
        }
        return outgoingShape;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkConnectableShape basicGetOutgoingShape() {
        return outgoingShape;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setOutgoingShape(ElkConnectableShape newOutgoingShape) {
        ElkConnectableShape oldOutgoingShape = outgoingShape;
        outgoingShape = newOutgoingShape;
        if (eNotificationRequired())
            eNotify(new ENotificationImpl(this, Notification.SET, ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SHAPE, oldOutgoingShape, outgoingShape));
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkConnectableShape getIncomingShape() {
        if (incomingShape != null && incomingShape.eIsProxy()) {
            InternalEObject oldIncomingShape = (InternalEObject)incomingShape;
            incomingShape = (ElkConnectableShape)eResolveProxy(oldIncomingShape);
            if (incomingShape != oldIncomingShape) {
                if (eNotificationRequired())
                    eNotify(new ENotificationImpl(this, Notification.RESOLVE, ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SHAPE, oldIncomingShape, incomingShape));
            }
        }
        return incomingShape;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkConnectableShape basicGetIncomingShape() {
        return incomingShape;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setIncomingShape(ElkConnectableShape newIncomingShape) {
        ElkConnectableShape oldIncomingShape = incomingShape;
        incomingShape = newIncomingShape;
        if (eNotificationRequired())
            eNotify(new ENotificationImpl(this, Notification.SET, ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SHAPE, oldIncomingShape, incomingShape));
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkEdgeSection> getOutgoingSections() {
        if (outgoingSections == null) {
            outgoingSections = new EObjectWithInverseResolvingEList.ManyInverse<ElkEdgeSection>(ElkEdgeSection.class, this, ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SECTIONS, ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SECTIONS);
        }
        return outgoingSections;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EList<ElkEdgeSection> getIncomingSections() {
        if (incomingSections == null) {
            incomingSections = new EObjectWithInverseResolvingEList.ManyInverse<ElkEdgeSection>(ElkEdgeSection.class, this, ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SECTIONS, ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SECTIONS);
        }
        return incomingSections;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public String getIdentifier() {
        return identifier;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setIdentifier(String newIdentifier) {
        String oldIdentifier = identifier;
        identifier = newIdentifier;
        if (eNotificationRequired())
            eNotify(new ENotificationImpl(this, Notification.SET, ElkGraphPackage.ELK_EDGE_SECTION__IDENTIFIER, oldIdentifier, identifier));
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setStartLocation(final double x, final double y) {
        setStartX(x);
        setStartY(y);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void setEndLocation(final double x, final double y) {
        setEndX(x);
        setEndY(y);
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
            case ElkGraphPackage.ELK_EDGE_SECTION__PARENT:
                if (eInternalContainer() != null)
                    msgs = eBasicRemoveFromContainer(msgs);
                return basicSetParent((ElkEdge)otherEnd, msgs);
            case ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SECTIONS:
                return ((InternalEList<InternalEObject>)(InternalEList<?>)getOutgoingSections()).basicAdd(otherEnd, msgs);
            case ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SECTIONS:
                return ((InternalEList<InternalEObject>)(InternalEList<?>)getIncomingSections()).basicAdd(otherEnd, msgs);
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
            case ElkGraphPackage.ELK_EDGE_SECTION__BEND_POINTS:
                return ((InternalEList<?>)getBendPoints()).basicRemove(otherEnd, msgs);
            case ElkGraphPackage.ELK_EDGE_SECTION__PARENT:
                return basicSetParent(null, msgs);
            case ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SECTIONS:
                return ((InternalEList<?>)getOutgoingSections()).basicRemove(otherEnd, msgs);
            case ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SECTIONS:
                return ((InternalEList<?>)getIncomingSections()).basicRemove(otherEnd, msgs);
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
            case ElkGraphPackage.ELK_EDGE_SECTION__PARENT:
                return eInternalContainer().eInverseRemove(this, ElkGraphPackage.ELK_EDGE__SECTIONS, ElkEdge.class, msgs);
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
            case ElkGraphPackage.ELK_EDGE_SECTION__START_X:
                return getStartX();
            case ElkGraphPackage.ELK_EDGE_SECTION__START_Y:
                return getStartY();
            case ElkGraphPackage.ELK_EDGE_SECTION__END_X:
                return getEndX();
            case ElkGraphPackage.ELK_EDGE_SECTION__END_Y:
                return getEndY();
            case ElkGraphPackage.ELK_EDGE_SECTION__BEND_POINTS:
                return getBendPoints();
            case ElkGraphPackage.ELK_EDGE_SECTION__PARENT:
                return getParent();
            case ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SHAPE:
                if (resolve) return getOutgoingShape();
                return basicGetOutgoingShape();
            case ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SHAPE:
                if (resolve) return getIncomingShape();
                return basicGetIncomingShape();
            case ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SECTIONS:
                return getOutgoingSections();
            case ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SECTIONS:
                return getIncomingSections();
            case ElkGraphPackage.ELK_EDGE_SECTION__IDENTIFIER:
                return getIdentifier();
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
            case ElkGraphPackage.ELK_EDGE_SECTION__START_X:
                setStartX((Double)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__START_Y:
                setStartY((Double)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__END_X:
                setEndX((Double)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__END_Y:
                setEndY((Double)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__BEND_POINTS:
                getBendPoints().clear();
                getBendPoints().addAll((Collection<? extends ElkBendPoint>)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__PARENT:
                setParent((ElkEdge)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SHAPE:
                setOutgoingShape((ElkConnectableShape)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SHAPE:
                setIncomingShape((ElkConnectableShape)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SECTIONS:
                getOutgoingSections().clear();
                getOutgoingSections().addAll((Collection<? extends ElkEdgeSection>)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SECTIONS:
                getIncomingSections().clear();
                getIncomingSections().addAll((Collection<? extends ElkEdgeSection>)newValue);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__IDENTIFIER:
                setIdentifier((String)newValue);
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
            case ElkGraphPackage.ELK_EDGE_SECTION__START_X:
                setStartX(START_X_EDEFAULT);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__START_Y:
                setStartY(START_Y_EDEFAULT);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__END_X:
                setEndX(END_X_EDEFAULT);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__END_Y:
                setEndY(END_Y_EDEFAULT);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__BEND_POINTS:
                getBendPoints().clear();
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__PARENT:
                setParent((ElkEdge)null);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SHAPE:
                setOutgoingShape((ElkConnectableShape)null);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SHAPE:
                setIncomingShape((ElkConnectableShape)null);
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SECTIONS:
                getOutgoingSections().clear();
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SECTIONS:
                getIncomingSections().clear();
                return;
            case ElkGraphPackage.ELK_EDGE_SECTION__IDENTIFIER:
                setIdentifier(IDENTIFIER_EDEFAULT);
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
            case ElkGraphPackage.ELK_EDGE_SECTION__START_X:
                return startX != START_X_EDEFAULT;
            case ElkGraphPackage.ELK_EDGE_SECTION__START_Y:
                return startY != START_Y_EDEFAULT;
            case ElkGraphPackage.ELK_EDGE_SECTION__END_X:
                return endX != END_X_EDEFAULT;
            case ElkGraphPackage.ELK_EDGE_SECTION__END_Y:
                return endY != END_Y_EDEFAULT;
            case ElkGraphPackage.ELK_EDGE_SECTION__BEND_POINTS:
                return bendPoints != null && !bendPoints.isEmpty();
            case ElkGraphPackage.ELK_EDGE_SECTION__PARENT:
                return getParent() != null;
            case ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SHAPE:
                return outgoingShape != null;
            case ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SHAPE:
                return incomingShape != null;
            case ElkGraphPackage.ELK_EDGE_SECTION__OUTGOING_SECTIONS:
                return outgoingSections != null && !outgoingSections.isEmpty();
            case ElkGraphPackage.ELK_EDGE_SECTION__INCOMING_SECTIONS:
                return incomingSections != null && !incomingSections.isEmpty();
            case ElkGraphPackage.ELK_EDGE_SECTION__IDENTIFIER:
                return IDENTIFIER_EDEFAULT == null ? identifier != null : !IDENTIFIER_EDEFAULT.equals(identifier);
        }
        return super.eIsSet(featureID);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public String toString() {
        if (eIsProxy()) return super.toString();

        StringBuffer result = new StringBuffer(super.toString());
        result.append(" (startX: ");
        result.append(startX);
        result.append(", startY: ");
        result.append(startY);
        result.append(", endX: ");
        result.append(endX);
        result.append(", endY: ");
        result.append(endY);
        result.append(", identifier: ");
        result.append(identifier);
        result.append(')');
        return result.toString();
    }

} //ElkEdgeSectionImpl
