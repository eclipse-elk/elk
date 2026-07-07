/**
 * ******************************************************************************
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
 *  ******************************************************************************
 */
package org.eclipse.elk.alg.graphviz.dot.dot.impl;

import org.eclipse.elk.alg.graphviz.dot.dot.DotPackage;
import org.eclipse.elk.alg.graphviz.dot.dot.EdgeOperator;
import org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget;
import org.eclipse.elk.alg.graphviz.dot.dot.Node;
import org.eclipse.elk.alg.graphviz.dot.dot.Subgraph;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;
import org.eclipse.emf.ecore.impl.MinimalEObjectImpl;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Edge Target</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.impl.EdgeTargetImpl#getOperator <em>Operator</em>}</li>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.impl.EdgeTargetImpl#getTargetSubgraph <em>Target Subgraph</em>}</li>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.impl.EdgeTargetImpl#getTargetnode <em>Targetnode</em>}</li>
 * </ul>
 *
 * @generated
 */
public class EdgeTargetImpl extends MinimalEObjectImpl.Container implements EdgeTarget
{
  /**
   * The default value of the '{@link #getOperator() <em>Operator</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getOperator()
   * @generated
   * @ordered
   */
  protected static final EdgeOperator OPERATOR_EDEFAULT = EdgeOperator.DIRECTED;

  /**
   * The cached value of the '{@link #getOperator() <em>Operator</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getOperator()
   * @generated
   * @ordered
   */
  protected EdgeOperator operator = OPERATOR_EDEFAULT;

  /**
   * The cached value of the '{@link #getTargetSubgraph() <em>Target Subgraph</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getTargetSubgraph()
   * @generated
   * @ordered
   */
  protected Subgraph targetSubgraph;

  /**
   * The cached value of the '{@link #getTargetnode() <em>Targetnode</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getTargetnode()
   * @generated
   * @ordered
   */
  protected Node targetnode;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected EdgeTargetImpl()
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
    return DotPackage.Literals.EDGE_TARGET;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EdgeOperator getOperator()
  {
    return operator;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setOperator(EdgeOperator newOperator)
  {
    EdgeOperator oldOperator = operator;
    operator = newOperator == null ? OPERATOR_EDEFAULT : newOperator;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, DotPackage.EDGE_TARGET__OPERATOR, oldOperator, operator));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public Subgraph getTargetSubgraph()
  {
    return targetSubgraph;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetTargetSubgraph(Subgraph newTargetSubgraph, NotificationChain msgs)
  {
    Subgraph oldTargetSubgraph = targetSubgraph;
    targetSubgraph = newTargetSubgraph;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, DotPackage.EDGE_TARGET__TARGET_SUBGRAPH, oldTargetSubgraph, newTargetSubgraph);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setTargetSubgraph(Subgraph newTargetSubgraph)
  {
    if (newTargetSubgraph != targetSubgraph)
    {
      NotificationChain msgs = null;
      if (targetSubgraph != null)
        msgs = ((InternalEObject)targetSubgraph).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - DotPackage.EDGE_TARGET__TARGET_SUBGRAPH, null, msgs);
      if (newTargetSubgraph != null)
        msgs = ((InternalEObject)newTargetSubgraph).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - DotPackage.EDGE_TARGET__TARGET_SUBGRAPH, null, msgs);
      msgs = basicSetTargetSubgraph(newTargetSubgraph, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, DotPackage.EDGE_TARGET__TARGET_SUBGRAPH, newTargetSubgraph, newTargetSubgraph));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public Node getTargetnode()
  {
    return targetnode;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetTargetnode(Node newTargetnode, NotificationChain msgs)
  {
    Node oldTargetnode = targetnode;
    targetnode = newTargetnode;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, DotPackage.EDGE_TARGET__TARGETNODE, oldTargetnode, newTargetnode);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setTargetnode(Node newTargetnode)
  {
    if (newTargetnode != targetnode)
    {
      NotificationChain msgs = null;
      if (targetnode != null)
        msgs = ((InternalEObject)targetnode).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - DotPackage.EDGE_TARGET__TARGETNODE, null, msgs);
      if (newTargetnode != null)
        msgs = ((InternalEObject)newTargetnode).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - DotPackage.EDGE_TARGET__TARGETNODE, null, msgs);
      msgs = basicSetTargetnode(newTargetnode, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, DotPackage.EDGE_TARGET__TARGETNODE, newTargetnode, newTargetnode));
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
      case DotPackage.EDGE_TARGET__TARGET_SUBGRAPH:
        return basicSetTargetSubgraph(null, msgs);
      case DotPackage.EDGE_TARGET__TARGETNODE:
        return basicSetTargetnode(null, msgs);
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
      case DotPackage.EDGE_TARGET__OPERATOR:
        return getOperator();
      case DotPackage.EDGE_TARGET__TARGET_SUBGRAPH:
        return getTargetSubgraph();
      case DotPackage.EDGE_TARGET__TARGETNODE:
        return getTargetnode();
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
      case DotPackage.EDGE_TARGET__OPERATOR:
        setOperator((EdgeOperator)newValue);
        return;
      case DotPackage.EDGE_TARGET__TARGET_SUBGRAPH:
        setTargetSubgraph((Subgraph)newValue);
        return;
      case DotPackage.EDGE_TARGET__TARGETNODE:
        setTargetnode((Node)newValue);
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
      case DotPackage.EDGE_TARGET__OPERATOR:
        setOperator(OPERATOR_EDEFAULT);
        return;
      case DotPackage.EDGE_TARGET__TARGET_SUBGRAPH:
        setTargetSubgraph((Subgraph)null);
        return;
      case DotPackage.EDGE_TARGET__TARGETNODE:
        setTargetnode((Node)null);
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
      case DotPackage.EDGE_TARGET__OPERATOR:
        return operator != OPERATOR_EDEFAULT;
      case DotPackage.EDGE_TARGET__TARGET_SUBGRAPH:
        return targetSubgraph != null;
      case DotPackage.EDGE_TARGET__TARGETNODE:
        return targetnode != null;
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
    result.append(" (operator: ");
    result.append(operator);
    result.append(')');
    return result.toString();
  }

} //EdgeTargetImpl
