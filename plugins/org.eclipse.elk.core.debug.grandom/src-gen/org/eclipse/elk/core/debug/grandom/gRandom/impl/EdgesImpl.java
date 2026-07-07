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
import org.eclipse.elk.core.debug.grandom.gRandom.Edges;
import org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;
import org.eclipse.emf.ecore.impl.MinimalEObjectImpl;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Edges</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl#isDensity <em>Density</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl#isTotal <em>Total</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl#isRelative <em>Relative</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl#isOutbound <em>Outbound</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl#getNEdges <em>NEdges</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl#isLabels <em>Labels</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.EdgesImpl#isSelfLoops <em>Self Loops</em>}</li>
 * </ul>
 *
 * @generated
 */
public class EdgesImpl extends MinimalEObjectImpl.Container implements Edges
{
  /**
   * The default value of the '{@link #isDensity() <em>Density</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isDensity()
   * @generated
   * @ordered
   */
  protected static final boolean DENSITY_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isDensity() <em>Density</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isDensity()
   * @generated
   * @ordered
   */
  protected boolean density = DENSITY_EDEFAULT;

  /**
   * The default value of the '{@link #isTotal() <em>Total</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isTotal()
   * @generated
   * @ordered
   */
  protected static final boolean TOTAL_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isTotal() <em>Total</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isTotal()
   * @generated
   * @ordered
   */
  protected boolean total = TOTAL_EDEFAULT;

  /**
   * The default value of the '{@link #isRelative() <em>Relative</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isRelative()
   * @generated
   * @ordered
   */
  protected static final boolean RELATIVE_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isRelative() <em>Relative</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isRelative()
   * @generated
   * @ordered
   */
  protected boolean relative = RELATIVE_EDEFAULT;

  /**
   * The default value of the '{@link #isOutbound() <em>Outbound</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isOutbound()
   * @generated
   * @ordered
   */
  protected static final boolean OUTBOUND_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isOutbound() <em>Outbound</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isOutbound()
   * @generated
   * @ordered
   */
  protected boolean outbound = OUTBOUND_EDEFAULT;

  /**
   * The cached value of the '{@link #getNEdges() <em>NEdges</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getNEdges()
   * @generated
   * @ordered
   */
  protected DoubleQuantity nEdges;

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
   * The default value of the '{@link #isSelfLoops() <em>Self Loops</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isSelfLoops()
   * @generated
   * @ordered
   */
  protected static final boolean SELF_LOOPS_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isSelfLoops() <em>Self Loops</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isSelfLoops()
   * @generated
   * @ordered
   */
  protected boolean selfLoops = SELF_LOOPS_EDEFAULT;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected EdgesImpl()
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
    return GRandomPackage.Literals.EDGES;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isDensity()
  {
    return density;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setDensity(boolean newDensity)
  {
    boolean oldDensity = density;
    density = newDensity;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.EDGES__DENSITY, oldDensity, density));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isTotal()
  {
    return total;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setTotal(boolean newTotal)
  {
    boolean oldTotal = total;
    total = newTotal;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.EDGES__TOTAL, oldTotal, total));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isRelative()
  {
    return relative;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setRelative(boolean newRelative)
  {
    boolean oldRelative = relative;
    relative = newRelative;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.EDGES__RELATIVE, oldRelative, relative));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isOutbound()
  {
    return outbound;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setOutbound(boolean newOutbound)
  {
    boolean oldOutbound = outbound;
    outbound = newOutbound;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.EDGES__OUTBOUND, oldOutbound, outbound));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public DoubleQuantity getNEdges()
  {
    return nEdges;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetNEdges(DoubleQuantity newNEdges, NotificationChain msgs)
  {
    DoubleQuantity oldNEdges = nEdges;
    nEdges = newNEdges;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, GRandomPackage.EDGES__NEDGES, oldNEdges, newNEdges);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setNEdges(DoubleQuantity newNEdges)
  {
    if (newNEdges != nEdges)
    {
      NotificationChain msgs = null;
      if (nEdges != null)
        msgs = ((InternalEObject)nEdges).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.EDGES__NEDGES, null, msgs);
      if (newNEdges != null)
        msgs = ((InternalEObject)newNEdges).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - GRandomPackage.EDGES__NEDGES, null, msgs);
      msgs = basicSetNEdges(newNEdges, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.EDGES__NEDGES, newNEdges, newNEdges));
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
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.EDGES__LABELS, oldLabels, labels));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isSelfLoops()
  {
    return selfLoops;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setSelfLoops(boolean newSelfLoops)
  {
    boolean oldSelfLoops = selfLoops;
    selfLoops = newSelfLoops;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.EDGES__SELF_LOOPS, oldSelfLoops, selfLoops));
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
      case GRandomPackage.EDGES__NEDGES:
        return basicSetNEdges(null, msgs);
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
      case GRandomPackage.EDGES__DENSITY:
        return isDensity();
      case GRandomPackage.EDGES__TOTAL:
        return isTotal();
      case GRandomPackage.EDGES__RELATIVE:
        return isRelative();
      case GRandomPackage.EDGES__OUTBOUND:
        return isOutbound();
      case GRandomPackage.EDGES__NEDGES:
        return getNEdges();
      case GRandomPackage.EDGES__LABELS:
        return isLabels();
      case GRandomPackage.EDGES__SELF_LOOPS:
        return isSelfLoops();
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
      case GRandomPackage.EDGES__DENSITY:
        setDensity((Boolean)newValue);
        return;
      case GRandomPackage.EDGES__TOTAL:
        setTotal((Boolean)newValue);
        return;
      case GRandomPackage.EDGES__RELATIVE:
        setRelative((Boolean)newValue);
        return;
      case GRandomPackage.EDGES__OUTBOUND:
        setOutbound((Boolean)newValue);
        return;
      case GRandomPackage.EDGES__NEDGES:
        setNEdges((DoubleQuantity)newValue);
        return;
      case GRandomPackage.EDGES__LABELS:
        setLabels((Boolean)newValue);
        return;
      case GRandomPackage.EDGES__SELF_LOOPS:
        setSelfLoops((Boolean)newValue);
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
      case GRandomPackage.EDGES__DENSITY:
        setDensity(DENSITY_EDEFAULT);
        return;
      case GRandomPackage.EDGES__TOTAL:
        setTotal(TOTAL_EDEFAULT);
        return;
      case GRandomPackage.EDGES__RELATIVE:
        setRelative(RELATIVE_EDEFAULT);
        return;
      case GRandomPackage.EDGES__OUTBOUND:
        setOutbound(OUTBOUND_EDEFAULT);
        return;
      case GRandomPackage.EDGES__NEDGES:
        setNEdges((DoubleQuantity)null);
        return;
      case GRandomPackage.EDGES__LABELS:
        setLabels(LABELS_EDEFAULT);
        return;
      case GRandomPackage.EDGES__SELF_LOOPS:
        setSelfLoops(SELF_LOOPS_EDEFAULT);
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
      case GRandomPackage.EDGES__DENSITY:
        return density != DENSITY_EDEFAULT;
      case GRandomPackage.EDGES__TOTAL:
        return total != TOTAL_EDEFAULT;
      case GRandomPackage.EDGES__RELATIVE:
        return relative != RELATIVE_EDEFAULT;
      case GRandomPackage.EDGES__OUTBOUND:
        return outbound != OUTBOUND_EDEFAULT;
      case GRandomPackage.EDGES__NEDGES:
        return nEdges != null;
      case GRandomPackage.EDGES__LABELS:
        return labels != LABELS_EDEFAULT;
      case GRandomPackage.EDGES__SELF_LOOPS:
        return selfLoops != SELF_LOOPS_EDEFAULT;
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
    result.append(" (density: ");
    result.append(density);
    result.append(", total: ");
    result.append(total);
    result.append(", relative: ");
    result.append(relative);
    result.append(", outbound: ");
    result.append(outbound);
    result.append(", labels: ");
    result.append(labels);
    result.append(", selfLoops: ");
    result.append(selfLoops);
    result.append(')');
    return result.toString();
  }

} //EdgesImpl
