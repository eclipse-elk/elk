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
import org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage;
import org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;
import org.eclipse.emf.ecore.impl.MinimalEObjectImpl;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Hierarchy</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.HierarchyImpl#getLevels <em>Levels</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.HierarchyImpl#getEdges <em>Edges</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.HierarchyImpl#getNumHierarchNodes <em>Num Hierarch Nodes</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.HierarchyImpl#getCrossHierarchRel <em>Cross Hierarch Rel</em>}</li>
 * </ul>
 *
 * @generated
 */
public class HierarchyImpl extends MinimalEObjectImpl.Container implements Hierarchy
{
  /**
   * The cached value of the '{@link #getLevels() <em>Levels</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getLevels()
   * @generated
   * @ordered
   */
  protected DoubleQuantity levels;

  /**
   * The cached value of the '{@link #getEdges() <em>Edges</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getEdges()
   * @generated
   * @ordered
   */
  protected DoubleQuantity edges;

  /**
   * The cached value of the '{@link #getNumHierarchNodes() <em>Num Hierarch Nodes</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getNumHierarchNodes()
   * @generated
   * @ordered
   */
  protected DoubleQuantity numHierarchNodes;

  /**
   * The cached value of the '{@link #getCrossHierarchRel() <em>Cross Hierarch Rel</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getCrossHierarchRel()
   * @generated
   * @ordered
   */
  protected DoubleQuantity crossHierarchRel;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected HierarchyImpl()
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
    return GRandomPackage.Literals.HIERARCHY;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public DoubleQuantity getLevels()
  {
    return levels;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetLevels(DoubleQuantity newLevels, NotificationChain msgs)
  {
    DoubleQuantity oldLevels = levels;
    levels = newLevels;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.HIERARCHY__LEVELS, oldLevels, newLevels);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setLevels(DoubleQuantity newLevels)
  {
    if (newLevels != levels)
    {
      NotificationChain msgs = null;
      if (levels != null)
        msgs = ((InternalEObject)levels).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.HIERARCHY__LEVELS, null, msgs);
      if (newLevels != null)
        msgs = ((InternalEObject)newLevels).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.HIERARCHY__LEVELS, null, msgs);
      msgs = basicSetLevels(newLevels, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.HIERARCHY__LEVELS, newLevels, newLevels));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public DoubleQuantity getEdges()
  {
    return edges;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetEdges(DoubleQuantity newEdges, NotificationChain msgs)
  {
    DoubleQuantity oldEdges = edges;
    edges = newEdges;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.HIERARCHY__EDGES, oldEdges, newEdges);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setEdges(DoubleQuantity newEdges)
  {
    if (newEdges != edges)
    {
      NotificationChain msgs = null;
      if (edges != null)
        msgs = ((InternalEObject)edges).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.HIERARCHY__EDGES, null, msgs);
      if (newEdges != null)
        msgs = ((InternalEObject)newEdges).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.HIERARCHY__EDGES, null, msgs);
      msgs = basicSetEdges(newEdges, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.HIERARCHY__EDGES, newEdges, newEdges));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public DoubleQuantity getNumHierarchNodes()
  {
    return numHierarchNodes;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetNumHierarchNodes(DoubleQuantity newNumHierarchNodes, NotificationChain msgs)
  {
    DoubleQuantity oldNumHierarchNodes = numHierarchNodes;
    numHierarchNodes = newNumHierarchNodes;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.HIERARCHY__NUM_HIERARCH_NODES, oldNumHierarchNodes, newNumHierarchNodes);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setNumHierarchNodes(DoubleQuantity newNumHierarchNodes)
  {
    if (newNumHierarchNodes != numHierarchNodes)
    {
      NotificationChain msgs = null;
      if (numHierarchNodes != null)
        msgs = ((InternalEObject)numHierarchNodes).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.HIERARCHY__NUM_HIERARCH_NODES, null, msgs);
      if (newNumHierarchNodes != null)
        msgs = ((InternalEObject)newNumHierarchNodes).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.HIERARCHY__NUM_HIERARCH_NODES, null, msgs);
      msgs = basicSetNumHierarchNodes(newNumHierarchNodes, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.HIERARCHY__NUM_HIERARCH_NODES, newNumHierarchNodes, newNumHierarchNodes));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public DoubleQuantity getCrossHierarchRel()
  {
    return crossHierarchRel;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetCrossHierarchRel(DoubleQuantity newCrossHierarchRel, NotificationChain msgs)
  {
    DoubleQuantity oldCrossHierarchRel = crossHierarchRel;
    crossHierarchRel = newCrossHierarchRel;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.HIERARCHY__CROSS_HIERARCH_REL, oldCrossHierarchRel, newCrossHierarchRel);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setCrossHierarchRel(DoubleQuantity newCrossHierarchRel)
  {
    if (newCrossHierarchRel != crossHierarchRel)
    {
      NotificationChain msgs = null;
      if (crossHierarchRel != null)
        msgs = ((InternalEObject)crossHierarchRel).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.HIERARCHY__CROSS_HIERARCH_REL, null, msgs);
      if (newCrossHierarchRel != null)
        msgs = ((InternalEObject)newCrossHierarchRel).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.HIERARCHY__CROSS_HIERARCH_REL, null, msgs);
      msgs = basicSetCrossHierarchRel(newCrossHierarchRel, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.HIERARCHY__CROSS_HIERARCH_REL, newCrossHierarchRel, newCrossHierarchRel));
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
      case GRandomPackage.HIERARCHY__LEVELS:
        return basicSetLevels(null, msgs);
      case GRandomPackage.HIERARCHY__EDGES:
        return basicSetEdges(null, msgs);
      case GRandomPackage.HIERARCHY__NUM_HIERARCH_NODES:
        return basicSetNumHierarchNodes(null, msgs);
      case GRandomPackage.HIERARCHY__CROSS_HIERARCH_REL:
        return basicSetCrossHierarchRel(null, msgs);
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
      case GRandomPackage.HIERARCHY__LEVELS:
        return getLevels();
      case GRandomPackage.HIERARCHY__EDGES:
        return getEdges();
      case GRandomPackage.HIERARCHY__NUM_HIERARCH_NODES:
        return getNumHierarchNodes();
      case GRandomPackage.HIERARCHY__CROSS_HIERARCH_REL:
        return getCrossHierarchRel();
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
      case GRandomPackage.HIERARCHY__LEVELS:
        setLevels((DoubleQuantity)newValue);
        return;
      case GRandomPackage.HIERARCHY__EDGES:
        setEdges((DoubleQuantity)newValue);
        return;
      case GRandomPackage.HIERARCHY__NUM_HIERARCH_NODES:
        setNumHierarchNodes((DoubleQuantity)newValue);
        return;
      case GRandomPackage.HIERARCHY__CROSS_HIERARCH_REL:
        setCrossHierarchRel((DoubleQuantity)newValue);
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
      case GRandomPackage.HIERARCHY__LEVELS:
        setLevels((DoubleQuantity)null);
        return;
      case GRandomPackage.HIERARCHY__EDGES:
        setEdges((DoubleQuantity)null);
        return;
      case GRandomPackage.HIERARCHY__NUM_HIERARCH_NODES:
        setNumHierarchNodes((DoubleQuantity)null);
        return;
      case GRandomPackage.HIERARCHY__CROSS_HIERARCH_REL:
        setCrossHierarchRel((DoubleQuantity)null);
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
      case GRandomPackage.HIERARCHY__LEVELS:
        return levels != null;
      case GRandomPackage.HIERARCHY__EDGES:
        return edges != null;
      case GRandomPackage.HIERARCHY__NUM_HIERARCH_NODES:
        return numHierarchNodes != null;
      case GRandomPackage.HIERARCHY__CROSS_HIERARCH_REL:
        return crossHierarchRel != null;
    }
    return super.eIsSet(featureID);
  }

} //HierarchyImpl
