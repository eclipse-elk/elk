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
package org.eclipse.elk.core.meta.metaData.impl;

import java.util.Collection;

import org.eclipse.elk.core.meta.metaData.MdAlgorithm;
import org.eclipse.elk.core.meta.metaData.MdCategory;
import org.eclipse.elk.core.meta.metaData.MdGraphFeature;
import org.eclipse.elk.core.meta.metaData.MdOptionSupport;
import org.eclipse.elk.core.meta.metaData.MetaDataPackage;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.common.util.EList;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;

import org.eclipse.emf.ecore.util.EDataTypeEList;
import org.eclipse.emf.ecore.util.EObjectContainmentEList;
import org.eclipse.emf.ecore.util.InternalEList;

import org.eclipse.xtext.common.types.JvmTypeReference;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Md Algorithm</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#isDeprecated <em>Deprecated</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#getProvider <em>Provider</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#getParameter <em>Parameter</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#getLabel <em>Label</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#getTargetClass <em>Target Class</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#getDescription <em>Description</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#getCategory <em>Category</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#getPreviewImage <em>Preview Image</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#getSupportedFeatures <em>Supported Features</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#getValidator <em>Validator</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl#getSupportedOptions <em>Supported Options</em>}</li>
 * </ul>
 *
 * @generated
 */
public class MdAlgorithmImpl extends MdBundleMemberImpl implements MdAlgorithm
{
  /**
   * The default value of the '{@link #isDeprecated() <em>Deprecated</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isDeprecated()
   * @generated
   * @ordered
   */
  protected static final boolean DEPRECATED_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isDeprecated() <em>Deprecated</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isDeprecated()
   * @generated
   * @ordered
   */
  protected boolean deprecated = DEPRECATED_EDEFAULT;

  /**
   * The cached value of the '{@link #getProvider() <em>Provider</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getProvider()
   * @generated
   * @ordered
   */
  protected JvmTypeReference provider;

  /**
   * The default value of the '{@link #getParameter() <em>Parameter</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getParameter()
   * @generated
   * @ordered
   */
  protected static final String PARAMETER_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getParameter() <em>Parameter</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getParameter()
   * @generated
   * @ordered
   */
  protected String parameter = PARAMETER_EDEFAULT;

  /**
   * The default value of the '{@link #getLabel() <em>Label</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getLabel()
   * @generated
   * @ordered
   */
  protected static final String LABEL_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getLabel() <em>Label</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getLabel()
   * @generated
   * @ordered
   */
  protected String label = LABEL_EDEFAULT;

  /**
   * The default value of the '{@link #getTargetClass() <em>Target Class</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getTargetClass()
   * @generated
   * @ordered
   */
  protected static final String TARGET_CLASS_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getTargetClass() <em>Target Class</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getTargetClass()
   * @generated
   * @ordered
   */
  protected String targetClass = TARGET_CLASS_EDEFAULT;

  /**
   * The default value of the '{@link #getDescription() <em>Description</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getDescription()
   * @generated
   * @ordered
   */
  protected static final String DESCRIPTION_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getDescription() <em>Description</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getDescription()
   * @generated
   * @ordered
   */
  protected String description = DESCRIPTION_EDEFAULT;

  /**
   * The cached value of the '{@link #getCategory() <em>Category</em>}' reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getCategory()
   * @generated
   * @ordered
   */
  protected MdCategory category;

  /**
   * The default value of the '{@link #getPreviewImage() <em>Preview Image</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getPreviewImage()
   * @generated
   * @ordered
   */
  protected static final String PREVIEW_IMAGE_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getPreviewImage() <em>Preview Image</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getPreviewImage()
   * @generated
   * @ordered
   */
  protected String previewImage = PREVIEW_IMAGE_EDEFAULT;

  /**
   * The cached value of the '{@link #getSupportedFeatures() <em>Supported Features</em>}' attribute list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getSupportedFeatures()
   * @generated
   * @ordered
   */
  protected EList<MdGraphFeature> supportedFeatures;

  /**
   * The cached value of the '{@link #getValidator() <em>Validator</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getValidator()
   * @generated
   * @ordered
   */
  protected JvmTypeReference validator;

  /**
   * The cached value of the '{@link #getSupportedOptions() <em>Supported Options</em>}' containment reference list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getSupportedOptions()
   * @generated
   * @ordered
   */
  protected EList<MdOptionSupport> supportedOptions;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected MdAlgorithmImpl()
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
    return MetaDataPackage.Literals.MD_ALGORITHM;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isDeprecated()
  {
    return deprecated;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setDeprecated(boolean newDeprecated)
  {
    boolean oldDeprecated = deprecated;
    deprecated = newDeprecated;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__DEPRECATED, oldDeprecated, deprecated));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public JvmTypeReference getProvider()
  {
    return provider;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetProvider(JvmTypeReference newProvider, NotificationChain msgs)
  {
    JvmTypeReference oldProvider = provider;
    provider = newProvider;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__PROVIDER, oldProvider, newProvider);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setProvider(JvmTypeReference newProvider)
  {
    if (newProvider != provider)
    {
      NotificationChain msgs = null;
      if (provider != null)
        msgs = ((InternalEObject)provider).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_ALGORITHM__PROVIDER, null, msgs);
      if (newProvider != null)
        msgs = ((InternalEObject)newProvider).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_ALGORITHM__PROVIDER, null, msgs);
      msgs = basicSetProvider(newProvider, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__PROVIDER, newProvider, newProvider));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getParameter()
  {
    return parameter;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setParameter(String newParameter)
  {
    String oldParameter = parameter;
    parameter = newParameter;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__PARAMETER, oldParameter, parameter));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getLabel()
  {
    return label;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setLabel(String newLabel)
  {
    String oldLabel = label;
    label = newLabel;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__LABEL, oldLabel, label));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getTargetClass()
  {
    return targetClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setTargetClass(String newTargetClass)
  {
    String oldTargetClass = targetClass;
    targetClass = newTargetClass;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__TARGET_CLASS, oldTargetClass, targetClass));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getDescription()
  {
    return description;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setDescription(String newDescription)
  {
    String oldDescription = description;
    description = newDescription;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__DESCRIPTION, oldDescription, description));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdCategory getCategory()
  {
    if (category != null && category.eIsProxy())
    {
      InternalEObject oldCategory = (InternalEObject)category;
      category = (MdCategory)eResolveProxy(oldCategory);
      if (category != oldCategory)
      {
        if (eNotificationRequired())
          eNotify(new ENotificationImpl(this, Notification.RESOLVE, MetaDataPackage.MD_ALGORITHM__CATEGORY, oldCategory, category));
      }
    }
    return category;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdCategory basicGetCategory()
  {
    return category;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setCategory(MdCategory newCategory)
  {
    MdCategory oldCategory = category;
    category = newCategory;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__CATEGORY, oldCategory, category));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getPreviewImage()
  {
    return previewImage;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setPreviewImage(String newPreviewImage)
  {
    String oldPreviewImage = previewImage;
    previewImage = newPreviewImage;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__PREVIEW_IMAGE, oldPreviewImage, previewImage));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EList<MdGraphFeature> getSupportedFeatures()
  {
    if (supportedFeatures == null)
    {
      supportedFeatures = new EDataTypeEList<MdGraphFeature>(MdGraphFeature.class, this, MetaDataPackage.MD_ALGORITHM__SUPPORTED_FEATURES);
    }
    return supportedFeatures;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public JvmTypeReference getValidator()
  {
    return validator;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetValidator(JvmTypeReference newValidator, NotificationChain msgs)
  {
    JvmTypeReference oldValidator = validator;
    validator = newValidator;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__VALIDATOR, oldValidator, newValidator);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setValidator(JvmTypeReference newValidator)
  {
    if (newValidator != validator)
    {
      NotificationChain msgs = null;
      if (validator != null)
        msgs = ((InternalEObject)validator).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_ALGORITHM__VALIDATOR, null, msgs);
      if (newValidator != null)
        msgs = ((InternalEObject)newValidator).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_ALGORITHM__VALIDATOR, null, msgs);
      msgs = basicSetValidator(newValidator, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_ALGORITHM__VALIDATOR, newValidator, newValidator));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EList<MdOptionSupport> getSupportedOptions()
  {
    if (supportedOptions == null)
    {
      supportedOptions = new EObjectContainmentEList<MdOptionSupport>(MdOptionSupport.class, this, MetaDataPackage.MD_ALGORITHM__SUPPORTED_OPTIONS);
    }
    return supportedOptions;
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
      case MetaDataPackage.MD_ALGORITHM__PROVIDER:
        return basicSetProvider(null, msgs);
      case MetaDataPackage.MD_ALGORITHM__VALIDATOR:
        return basicSetValidator(null, msgs);
      case MetaDataPackage.MD_ALGORITHM__SUPPORTED_OPTIONS:
        return ((InternalEList<?>)getSupportedOptions()).basicRemove(otherEnd, msgs);
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
      case MetaDataPackage.MD_ALGORITHM__DEPRECATED:
        return isDeprecated();
      case MetaDataPackage.MD_ALGORITHM__PROVIDER:
        return getProvider();
      case MetaDataPackage.MD_ALGORITHM__PARAMETER:
        return getParameter();
      case MetaDataPackage.MD_ALGORITHM__LABEL:
        return getLabel();
      case MetaDataPackage.MD_ALGORITHM__TARGET_CLASS:
        return getTargetClass();
      case MetaDataPackage.MD_ALGORITHM__DESCRIPTION:
        return getDescription();
      case MetaDataPackage.MD_ALGORITHM__CATEGORY:
        if (resolve) return getCategory();
        return basicGetCategory();
      case MetaDataPackage.MD_ALGORITHM__PREVIEW_IMAGE:
        return getPreviewImage();
      case MetaDataPackage.MD_ALGORITHM__SUPPORTED_FEATURES:
        return getSupportedFeatures();
      case MetaDataPackage.MD_ALGORITHM__VALIDATOR:
        return getValidator();
      case MetaDataPackage.MD_ALGORITHM__SUPPORTED_OPTIONS:
        return getSupportedOptions();
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
      case MetaDataPackage.MD_ALGORITHM__DEPRECATED:
        setDeprecated((Boolean)newValue);
        return;
      case MetaDataPackage.MD_ALGORITHM__PROVIDER:
        setProvider((JvmTypeReference)newValue);
        return;
      case MetaDataPackage.MD_ALGORITHM__PARAMETER:
        setParameter((String)newValue);
        return;
      case MetaDataPackage.MD_ALGORITHM__LABEL:
        setLabel((String)newValue);
        return;
      case MetaDataPackage.MD_ALGORITHM__TARGET_CLASS:
        setTargetClass((String)newValue);
        return;
      case MetaDataPackage.MD_ALGORITHM__DESCRIPTION:
        setDescription((String)newValue);
        return;
      case MetaDataPackage.MD_ALGORITHM__CATEGORY:
        setCategory((MdCategory)newValue);
        return;
      case MetaDataPackage.MD_ALGORITHM__PREVIEW_IMAGE:
        setPreviewImage((String)newValue);
        return;
      case MetaDataPackage.MD_ALGORITHM__SUPPORTED_FEATURES:
        getSupportedFeatures().clear();
        getSupportedFeatures().addAll((Collection<? extends MdGraphFeature>)newValue);
        return;
      case MetaDataPackage.MD_ALGORITHM__VALIDATOR:
        setValidator((JvmTypeReference)newValue);
        return;
      case MetaDataPackage.MD_ALGORITHM__SUPPORTED_OPTIONS:
        getSupportedOptions().clear();
        getSupportedOptions().addAll((Collection<? extends MdOptionSupport>)newValue);
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
      case MetaDataPackage.MD_ALGORITHM__DEPRECATED:
        setDeprecated(DEPRECATED_EDEFAULT);
        return;
      case MetaDataPackage.MD_ALGORITHM__PROVIDER:
        setProvider((JvmTypeReference)null);
        return;
      case MetaDataPackage.MD_ALGORITHM__PARAMETER:
        setParameter(PARAMETER_EDEFAULT);
        return;
      case MetaDataPackage.MD_ALGORITHM__LABEL:
        setLabel(LABEL_EDEFAULT);
        return;
      case MetaDataPackage.MD_ALGORITHM__TARGET_CLASS:
        setTargetClass(TARGET_CLASS_EDEFAULT);
        return;
      case MetaDataPackage.MD_ALGORITHM__DESCRIPTION:
        setDescription(DESCRIPTION_EDEFAULT);
        return;
      case MetaDataPackage.MD_ALGORITHM__CATEGORY:
        setCategory((MdCategory)null);
        return;
      case MetaDataPackage.MD_ALGORITHM__PREVIEW_IMAGE:
        setPreviewImage(PREVIEW_IMAGE_EDEFAULT);
        return;
      case MetaDataPackage.MD_ALGORITHM__SUPPORTED_FEATURES:
        getSupportedFeatures().clear();
        return;
      case MetaDataPackage.MD_ALGORITHM__VALIDATOR:
        setValidator((JvmTypeReference)null);
        return;
      case MetaDataPackage.MD_ALGORITHM__SUPPORTED_OPTIONS:
        getSupportedOptions().clear();
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
      case MetaDataPackage.MD_ALGORITHM__DEPRECATED:
        return deprecated != DEPRECATED_EDEFAULT;
      case MetaDataPackage.MD_ALGORITHM__PROVIDER:
        return provider != null;
      case MetaDataPackage.MD_ALGORITHM__PARAMETER:
        return PARAMETER_EDEFAULT == null ? parameter != null : !PARAMETER_EDEFAULT.equals(parameter);
      case MetaDataPackage.MD_ALGORITHM__LABEL:
        return LABEL_EDEFAULT == null ? label != null : !LABEL_EDEFAULT.equals(label);
      case MetaDataPackage.MD_ALGORITHM__TARGET_CLASS:
        return TARGET_CLASS_EDEFAULT == null ? targetClass != null : !TARGET_CLASS_EDEFAULT.equals(targetClass);
      case MetaDataPackage.MD_ALGORITHM__DESCRIPTION:
        return DESCRIPTION_EDEFAULT == null ? description != null : !DESCRIPTION_EDEFAULT.equals(description);
      case MetaDataPackage.MD_ALGORITHM__CATEGORY:
        return category != null;
      case MetaDataPackage.MD_ALGORITHM__PREVIEW_IMAGE:
        return PREVIEW_IMAGE_EDEFAULT == null ? previewImage != null : !PREVIEW_IMAGE_EDEFAULT.equals(previewImage);
      case MetaDataPackage.MD_ALGORITHM__SUPPORTED_FEATURES:
        return supportedFeatures != null && !supportedFeatures.isEmpty();
      case MetaDataPackage.MD_ALGORITHM__VALIDATOR:
        return validator != null;
      case MetaDataPackage.MD_ALGORITHM__SUPPORTED_OPTIONS:
        return supportedOptions != null && !supportedOptions.isEmpty();
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
    result.append(" (deprecated: ");
    result.append(deprecated);
    result.append(", parameter: ");
    result.append(parameter);
    result.append(", label: ");
    result.append(label);
    result.append(", targetClass: ");
    result.append(targetClass);
    result.append(", description: ");
    result.append(description);
    result.append(", previewImage: ");
    result.append(previewImage);
    result.append(", supportedFeatures: ");
    result.append(supportedFeatures);
    result.append(')');
    return result.toString();
  }

} //MdAlgorithmImpl
