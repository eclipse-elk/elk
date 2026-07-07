/**
 * ******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *  ******************************************************************************
 */
package org.eclipse.elk.core.debug.grandom.gRandom.impl;

import org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity;
import org.eclipse.elk.core.debug.grandom.gRandom.Flow;
import org.eclipse.elk.core.debug.grandom.gRandom.FlowType;
import org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage;
import org.eclipse.elk.core.debug.grandom.gRandom.Side;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;
import org.eclipse.emf.ecore.impl.MinimalEObjectImpl;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Flow</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.FlowImpl#getFlowType <em>Flow Type</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.FlowImpl#getSide <em>Side</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.FlowImpl#getAmount <em>Amount</em>}</li>
 * </ul>
 *
 * @generated
 */
public class FlowImpl extends MinimalEObjectImpl.Container implements Flow
{
  /**
   * The default value of the '{@link #getFlowType() <em>Flow Type</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getFlowType()
   * @generated
   * @ordered
   */
  protected static final FlowType FLOW_TYPE_EDEFAULT = FlowType.INCOMING;

  /**
   * The cached value of the '{@link #getFlowType() <em>Flow Type</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getFlowType()
   * @generated
   * @ordered
   */
  protected FlowType flowType = FLOW_TYPE_EDEFAULT;

  /**
   * The default value of the '{@link #getSide() <em>Side</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getSide()
   * @generated
   * @ordered
   */
  protected static final Side SIDE_EDEFAULT = Side.NORTH;

  /**
   * The cached value of the '{@link #getSide() <em>Side</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getSide()
   * @generated
   * @ordered
   */
  protected Side side = SIDE_EDEFAULT;

  /**
   * The cached value of the '{@link #getAmount() <em>Amount</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getAmount()
   * @generated
   * @ordered
   */
  protected DoubleQuantity amount;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected FlowImpl()
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
    return GRandomPackage.Literals.FLOW;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public FlowType getFlowType()
  {
    return flowType;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setFlowType(FlowType newFlowType)
  {
    FlowType oldFlowType = flowType;
    flowType = newFlowType == null ? FLOW_TYPE_EDEFAULT : newFlowType;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.FLOW__FLOW_TYPE, oldFlowType, flowType));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public Side getSide()
  {
    return side;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setSide(Side newSide)
  {
    Side oldSide = side;
    side = newSide == null ? SIDE_EDEFAULT : newSide;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.FLOW__SIDE, oldSide, side));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public DoubleQuantity getAmount()
  {
    return amount;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetAmount(DoubleQuantity newAmount, NotificationChain msgs)
  {
    DoubleQuantity oldAmount = amount;
    amount = newAmount;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.FLOW__AMOUNT, oldAmount, newAmount);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setAmount(DoubleQuantity newAmount)
  {
    if (newAmount != amount)
    {
      NotificationChain msgs = null;
      if (amount != null)
        msgs = ((InternalEObject)amount).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.FLOW__AMOUNT, null, msgs);
      if (newAmount != null)
        msgs = ((InternalEObject)newAmount).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.FLOW__AMOUNT, null, msgs);
      msgs = basicSetAmount(newAmount, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.FLOW__AMOUNT, newAmount, newAmount));
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
      case GRandomPackage.FLOW__AMOUNT:
        return basicSetAmount(null, msgs);
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
      case GRandomPackage.FLOW__FLOW_TYPE:
        return getFlowType();
      case GRandomPackage.FLOW__SIDE:
        return getSide();
      case GRandomPackage.FLOW__AMOUNT:
        return getAmount();
    }
    return super.eGet(featureID, resolve, coreType);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public void eSet(int featureID, Object newValue)
  {
    switch (featureID)
    {
      case GRandomPackage.FLOW__FLOW_TYPE:
        setFlowType((FlowType)newValue);
        return;
      case GRandomPackage.FLOW__SIDE:
        setSide((Side)newValue);
        return;
      case GRandomPackage.FLOW__AMOUNT:
        setAmount((DoubleQuantity)newValue);
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
      case GRandomPackage.FLOW__FLOW_TYPE:
        setFlowType(FLOW_TYPE_EDEFAULT);
        return;
      case GRandomPackage.FLOW__SIDE:
        setSide(SIDE_EDEFAULT);
        return;
      case GRandomPackage.FLOW__AMOUNT:
        setAmount((DoubleQuantity)null);
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
      case GRandomPackage.FLOW__FLOW_TYPE:
        return flowType != FLOW_TYPE_EDEFAULT;
      case GRandomPackage.FLOW__SIDE:
        return side != SIDE_EDEFAULT;
      case GRandomPackage.FLOW__AMOUNT:
        return amount != null;
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
    result.append(" (flowType: ");
    result.append(flowType);
    result.append(", side: ");
    result.append(side);
    result.append(')');
    return result.toString();
  }

} //FlowImpl
