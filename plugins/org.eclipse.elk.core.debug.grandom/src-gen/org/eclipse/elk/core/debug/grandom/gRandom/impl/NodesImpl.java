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

import org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity;
import org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage;
import org.eclipse.elk.core.debug.grandom.gRandom.Nodes;
import org.eclipse.elk.core.debug.grandom.gRandom.Ports;
import org.eclipse.elk.core.debug.grandom.gRandom.Size;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;
import org.eclipse.emf.ecore.impl.MinimalEObjectImpl;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Nodes</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.NodesImpl#getNNodes <em>NNodes</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.NodesImpl#getPorts <em>Ports</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.NodesImpl#isLabels <em>Labels</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.NodesImpl#getSize <em>Size</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.NodesImpl#isRemoveIsolated <em>Remove Isolated</em>}</li>
 * </ul>
 *
 * @generated
 */
public class NodesImpl extends MinimalEObjectImpl.Container implements Nodes
{
  /**
   * The cached value of the '{@link #getNNodes() <em>NNodes</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getNNodes()
   * @generated
   * @ordered
   */
  protected DoubleQuantity nNodes;

  /**
   * The cached value of the '{@link #getPorts() <em>Ports</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getPorts()
   * @generated
   * @ordered
   */
  protected Ports ports;

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
   * The cached value of the '{@link #getSize() <em>Size</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getSize()
   * @generated
   * @ordered
   */
  protected Size size;

  /**
   * The default value of the '{@link #isRemoveIsolated() <em>Remove Isolated</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isRemoveIsolated()
   * @generated
   * @ordered
   */
  protected static final boolean REMOVE_ISOLATED_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isRemoveIsolated() <em>Remove Isolated</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isRemoveIsolated()
   * @generated
   * @ordered
   */
  protected boolean removeIsolated = REMOVE_ISOLATED_EDEFAULT;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected NodesImpl()
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
    return GRandomPackage.Literals.NODES;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public DoubleQuantity getNNodes()
  {
    return nNodes;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetNNodes(DoubleQuantity newNNodes, NotificationChain msgs)
  {
    DoubleQuantity oldNNodes = nNodes;
    nNodes = newNNodes;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.NODES__NNODES, oldNNodes, newNNodes);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setNNodes(DoubleQuantity newNNodes)
  {
    if (newNNodes != nNodes)
    {
      NotificationChain msgs = null;
      if (nNodes != null)
        msgs = ((InternalEObject)nNodes).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.NODES__NNODES, null, msgs);
      if (newNNodes != null)
        msgs = ((InternalEObject)newNNodes).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.NODES__NNODES, null, msgs);
      msgs = basicSetNNodes(newNNodes, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.NODES__NNODES, newNNodes, newNNodes));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public Ports getPorts()
  {
    return ports;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetPorts(Ports newPorts, NotificationChain msgs)
  {
    Ports oldPorts = ports;
    ports = newPorts;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.NODES__PORTS, oldPorts, newPorts);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setPorts(Ports newPorts)
  {
    if (newPorts != ports)
    {
      NotificationChain msgs = null;
      if (ports != null)
        msgs = ((InternalEObject)ports).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.NODES__PORTS, null, msgs);
      if (newPorts != null)
        msgs = ((InternalEObject)newPorts).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.NODES__PORTS, null, msgs);
      msgs = basicSetPorts(newPorts, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.NODES__PORTS, newPorts, newPorts));
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
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.NODES__LABELS, oldLabels, labels));
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
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.NODES__SIZE, oldSize, newSize);
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
        msgs = ((InternalEObject)size).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.NODES__SIZE, null, msgs);
      if (newSize != null)
        msgs = ((InternalEObject)newSize).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.NODES__SIZE, null, msgs);
      msgs = basicSetSize(newSize, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.NODES__SIZE, newSize, newSize));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isRemoveIsolated()
  {
    return removeIsolated;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setRemoveIsolated(boolean newRemoveIsolated)
  {
    boolean oldRemoveIsolated = removeIsolated;
    removeIsolated = newRemoveIsolated;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.NODES__REMOVE_ISOLATED, oldRemoveIsolated, removeIsolated));
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
      case GRandomPackage.NODES__NNODES:
        return basicSetNNodes(null, msgs);
      case GRandomPackage.NODES__PORTS:
        return basicSetPorts(null, msgs);
      case GRandomPackage.NODES__SIZE:
        return basicSetSize(null, msgs);
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
      case GRandomPackage.NODES__NNODES:
        return getNNodes();
      case GRandomPackage.NODES__PORTS:
        return getPorts();
      case GRandomPackage.NODES__LABELS:
        return isLabels();
      case GRandomPackage.NODES__SIZE:
        return getSize();
      case GRandomPackage.NODES__REMOVE_ISOLATED:
        return isRemoveIsolated();
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
      case GRandomPackage.NODES__NNODES:
        setNNodes((DoubleQuantity)newValue);
        return;
      case GRandomPackage.NODES__PORTS:
        setPorts((Ports)newValue);
        return;
      case GRandomPackage.NODES__LABELS:
        setLabels((Boolean)newValue);
        return;
      case GRandomPackage.NODES__SIZE:
        setSize((Size)newValue);
        return;
      case GRandomPackage.NODES__REMOVE_ISOLATED:
        setRemoveIsolated((Boolean)newValue);
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
      case GRandomPackage.NODES__NNODES:
        setNNodes((DoubleQuantity)null);
        return;
      case GRandomPackage.NODES__PORTS:
        setPorts((Ports)null);
        return;
      case GRandomPackage.NODES__LABELS:
        setLabels(LABELS_EDEFAULT);
        return;
      case GRandomPackage.NODES__SIZE:
        setSize((Size)null);
        return;
      case GRandomPackage.NODES__REMOVE_ISOLATED:
        setRemoveIsolated(REMOVE_ISOLATED_EDEFAULT);
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
      case GRandomPackage.NODES__NNODES:
        return nNodes != null;
      case GRandomPackage.NODES__PORTS:
        return ports != null;
      case GRandomPackage.NODES__LABELS:
        return labels != LABELS_EDEFAULT;
      case GRandomPackage.NODES__SIZE:
        return size != null;
      case GRandomPackage.NODES__REMOVE_ISOLATED:
        return removeIsolated != REMOVE_ISOLATED_EDEFAULT;
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
    result.append(", removeIsolated: ");
    result.append(removeIsolated);
    result.append(')');
    return result.toString();
  }

} //NodesImpl
