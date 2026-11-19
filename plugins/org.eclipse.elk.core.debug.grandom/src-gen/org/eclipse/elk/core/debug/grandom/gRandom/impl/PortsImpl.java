/**
 * ******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
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
 *  ******************************************************************************
 */
package org.eclipse.elk.core.debug.grandom.gRandom.impl;

import java.util.Collection;

import org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType;
import org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity;
import org.eclipse.elk.core.debug.grandom.gRandom.Flow;
import org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage;
import org.eclipse.elk.core.debug.grandom.gRandom.Ports;
import org.eclipse.elk.core.debug.grandom.gRandom.Size;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.common.util.EList;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;
import org.eclipse.emf.ecore.impl.MinimalEObjectImpl;

import org.eclipse.emf.ecore.util.EObjectContainmentEList;
import org.eclipse.emf.ecore.util.InternalEList;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Ports</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.PortsImpl#isLabels <em>Labels</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.PortsImpl#getReUse <em>Re Use</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.PortsImpl#getSize <em>Size</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.PortsImpl#getConstraint <em>Constraint</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.PortsImpl#getFlow <em>Flow</em>}</li>
 * </ul>
 *
 * @generated
 */
public class PortsImpl extends MinimalEObjectImpl.Container implements Ports
{
  /**
   * The default value of the '{@link #isLabels() <em>Labels</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isLabels()
   * @generated
   * @ordered
   */
  protected static final boolean LABELS_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isLabels() <em>Labels</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isLabels()
   * @generated
   * @ordered
   */
  protected boolean labels = LABELS_EDEFAULT;

  /**
   * The cached value of the '{@link #getReUse() <em>Re Use</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getReUse()
   * @generated
   * @ordered
   */
  protected DoubleQuantity reUse;

  /**
   * The cached value of the '{@link #getSize() <em>Size</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getSize()
   * @generated
   * @ordered
   */
  protected Size size;

  /**
   * The default value of the '{@link #getConstraint() <em>Constraint</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getConstraint()
   * @generated
   * @ordered
   */
  protected static final ConstraintType CONSTRAINT_EDEFAULT = ConstraintType.FREE;

  /**
   * The cached value of the '{@link #getConstraint() <em>Constraint</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getConstraint()
   * @generated
   * @ordered
   */
  protected ConstraintType constraint = CONSTRAINT_EDEFAULT;

  /**
   * The cached value of the '{@link #getFlow() <em>Flow</em>}' containment reference list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getFlow()
   * @generated
   * @ordered
   */
  protected EList<Flow> flow;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected PortsImpl()
  {
    super();
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  protected EClass eStaticClass()
  {
    return GRandomPackage.Literals.PORTS;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isLabels()
  {
    return labels;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setLabels(boolean newLabels)
  {
    boolean oldLabels = labels;
    labels = newLabels;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.PORTS__LABELS, oldLabels, labels));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public DoubleQuantity getReUse()
  {
    return reUse;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetReUse(DoubleQuantity newReUse, NotificationChain msgs)
  {
    DoubleQuantity oldReUse = reUse;
    reUse = newReUse;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.PORTS__RE_USE, oldReUse, newReUse);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setReUse(DoubleQuantity newReUse)
  {
    if (newReUse != reUse)
    {
      NotificationChain msgs = null;
      if (reUse != null)
        msgs = ((InternalEObject)reUse).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.PORTS__RE_USE, null, msgs);
      if (newReUse != null)
        msgs = ((InternalEObject)newReUse).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.PORTS__RE_USE, null, msgs);
      msgs = basicSetReUse(newReUse, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.PORTS__RE_USE, newReUse, newReUse));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public Size getSize()
  {
    return size;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetSize(Size newSize, NotificationChain msgs)
  {
    Size oldSize = size;
    size = newSize;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.PORTS__SIZE, oldSize, newSize);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setSize(Size newSize)
  {
    if (newSize != size)
    {
      NotificationChain msgs = null;
      if (size != null)
        msgs = ((InternalEObject)size).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.PORTS__SIZE, null, msgs);
      if (newSize != null)
        msgs = ((InternalEObject)newSize).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.PORTS__SIZE, null, msgs);
      msgs = basicSetSize(newSize, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.PORTS__SIZE, newSize, newSize));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public ConstraintType getConstraint()
  {
    return constraint;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setConstraint(ConstraintType newConstraint)
  {
    ConstraintType oldConstraint = constraint;
    constraint = newConstraint == null ? CONSTRAINT_EDEFAULT : newConstraint;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.PORTS__CONSTRAINT, oldConstraint, constraint));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EList<Flow> getFlow()
  {
    if (flow == null)
    {
      flow = new EObjectContainmentEList<Flow>(Flow.class, this, GRandomPackage.PORTS__FLOW);
    }
    return flow;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public NotificationChain eInverseRemove(InternalEObject otherEnd, int featureID, NotificationChain msgs)
  {
    switch (featureID)
    {
      case GRandomPackage.PORTS__RE_USE:
        return basicSetReUse(null, msgs);
      case GRandomPackage.PORTS__SIZE:
        return basicSetSize(null, msgs);
      case GRandomPackage.PORTS__FLOW:
        return ((InternalEList<?>)getFlow()).basicRemove(otherEnd, msgs);
    }
    return super.eInverseRemove(otherEnd, featureID, msgs);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public Object eGet(int featureID, boolean resolve, boolean coreType)
  {
    switch (featureID)
    {
      case GRandomPackage.PORTS__LABELS:
        return isLabels();
      case GRandomPackage.PORTS__RE_USE:
        return getReUse();
      case GRandomPackage.PORTS__SIZE:
        return getSize();
      case GRandomPackage.PORTS__CONSTRAINT:
        return getConstraint();
      case GRandomPackage.PORTS__FLOW:
        return getFlow();
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
  public void eSet(int featureID, Object newValue)
  {
    switch (featureID)
    {
      case GRandomPackage.PORTS__LABELS:
        setLabels((Boolean)newValue);
        return;
      case GRandomPackage.PORTS__RE_USE:
        setReUse((DoubleQuantity)newValue);
        return;
      case GRandomPackage.PORTS__SIZE:
        setSize((Size)newValue);
        return;
      case GRandomPackage.PORTS__CONSTRAINT:
        setConstraint((ConstraintType)newValue);
        return;
      case GRandomPackage.PORTS__FLOW:
        getFlow().clear();
        getFlow().addAll((Collection<? extends Flow>)newValue);
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
  public void eUnset(int featureID)
  {
    switch (featureID)
    {
      case GRandomPackage.PORTS__LABELS:
        setLabels(LABELS_EDEFAULT);
        return;
      case GRandomPackage.PORTS__RE_USE:
        setReUse((DoubleQuantity)null);
        return;
      case GRandomPackage.PORTS__SIZE:
        setSize((Size)null);
        return;
      case GRandomPackage.PORTS__CONSTRAINT:
        setConstraint(CONSTRAINT_EDEFAULT);
        return;
      case GRandomPackage.PORTS__FLOW:
        getFlow().clear();
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
  public boolean eIsSet(int featureID)
  {
    switch (featureID)
    {
      case GRandomPackage.PORTS__LABELS:
        return labels != LABELS_EDEFAULT;
      case GRandomPackage.PORTS__RE_USE:
        return reUse != null;
      case GRandomPackage.PORTS__SIZE:
        return size != null;
      case GRandomPackage.PORTS__CONSTRAINT:
        return constraint != CONSTRAINT_EDEFAULT;
      case GRandomPackage.PORTS__FLOW:
        return flow != null && !flow.isEmpty();
    }
    return super.eIsSet(featureID);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public String toString()
  {
    if (eIsProxy()) return super.toString();

    StringBuilder result = new StringBuilder(super.toString());
    result.append(" (labels: ");
    result.append(labels);
    result.append(", constraint: ");
    result.append(constraint);
    result.append(')');
    return result.toString();
  }

} //PortsImpl
